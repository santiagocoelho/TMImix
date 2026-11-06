classdef TMImix
    % TMImix  Composable, baseline-normalized LTE diffusion MRI model.
    % Rows of the protocol are [b, Delta, delta] (ms/um^2, ms, ms).
    properties (SetAccess=private)
        protocol
        parts
        exchange
        steps
        mu
        weights
        fractions
    end
    methods
        function obj = TMImix(protocol,parts,opt)
            arguments
                protocol (:,3) double {mustBeReal,mustBeFinite}
                parts (1,:) string
                opt.Exchange (1,1) string {mustBeMember(opt.Exchange,["finite","narrow"])} = "finite"
                opt.Steps (1,1) double {mustBeInteger,mustBePositive} = 16
                opt.Order (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opt.Order,4)} = 48
            end
            keys = ["stick","ecs","permeableStick","sphere","dot","water"];
            names = ["fStick","fECS","fPerm","fSphere","fDot","fWater"];
            assert(~isempty(parts) && all(ismember(parts,keys)) && ...
                numel(unique(parts))==numel(parts),'Invalid or duplicate compartment.');
            assert(~ismember("permeableStick",parts) || ismember("ecs",parts), ...
                'permeableStick requires its exchange partner ecs.');
            assert(~isempty(protocol) && all(protocol(:,1)>=0 & protocol(:,2)>0 & ...
                protocol(:,3)>=0 & protocol(:,3)<=protocol(:,2)), ...
                'Expected b >= 0, Delta > 0 and 0 <= delta <= Delta.');
            obj.protocol=protocol; obj.parts=parts; obj.exchange=opt.Exchange;
            obj.steps=opt.Steps;
            [~,ix]=ismember(parts,keys); obj.fractions=names(ix);
            % Gauss-Legendre rule on [0,1].
            j=1:opt.Order-1; v=j./sqrt(4*j.^2-1);
            [V,L]=eig(diag(v,1)+diag(v,-1),'vector');
            [x,i]=sort(L); obj.mu=(x'+1)/2; obj.weights=V(1,i)'.^2;
        end

        function p = defaults(obj)
            % Scalar physical parameters; replace any scalar by [lo hi] in a fit spec.
            p=struct;
            for k=obj.fractions, p.(k)=1/numel(obj.parts); end
            if isfield(p,'fStick'), p.Dstick=2; end
            if isfield(p,'fECS'), p.DePar=1.5; p.DePerp=0.7; end
            if isfield(p,'fPerm'), p.Dperm=2; p.tau=30; end
            if isfield(p,'fSphere'), p.radius=7; p.Dsoma=2; end
            if isfield(p,'fWater'), p.Dwater=3; end
            if any(ismember(obj.parts,["stick","ecs","permeableStick"])), p.p2=0.6; end
        end

        function bounds = freeBounds(obj)
            % Editable simulation/fit ranges for every active parameter.
            ranges = struct('Dstick',[0.5 3], 'DePar',[0.3 2.5], 'DePerp',[0.1 1.5], ...
                'Dperm',[0.5 3], 'tau',[5 100], 'radius',[3 12], 'Dsoma',[0.5 3], ...
                'Dwater',[2 3.5], 'p2',[0 1]);
            bounds = obj.defaults();
            for name = string(fieldnames(bounds))'
                if ismember(name,obj.fractions)
                    bounds.(name) = [0 1];
                else
                    bounds.(name) = ranges.(name);
                end
            end
        end

        function K = kernel(obj,p,mu)
            % Full angular response, shells x direction cosines (no SH truncation).
            tmix.validate(obj,p,false);
            validateattributes(mu,{'numeric'},{'real','finite','vector','nonempty'});
            assert(all(abs(mu)<=1+1e-12),'Direction cosines must lie in [-1,1].');
            K=obj.response(p,min(1,abs(reshape(mu,1,[]))).^2);
        end

        function S = invariants(obj,p)
            % Columns: S_0 = K_0, S_2 = p2 * abs(K_2), normalized by b=0.
            K=obj.kernel(p,obj.mu);
            S=[K*obj.weights, zeros(size(K,1),1)];
            if isfield(p,'p2')
                S(:,2)=p.p2*abs(K*(obj.weights.*(3*obj.mu'.^2-1)/2));
            end
            S(obj.protocol(:,1)==0,:)=repmat([1 0],nnz(obj.protocol(:,1)==0),1);
        end

        function Y = signal(obj,p,dirs,axes,mass)
            % Each fODF node has probability mass, including quadrature area.
            % Returns directions x shells. p.p2 is ignored: the fODF sets it.
            tmix.validate(obj,p,false); tmix.unitRows(dirs);
            [axes,mass]=tmix.odf(axes,mass);
            Y=zeros(size(dirs,1),size(obj.protocol,1));
            % Block both axes and gradients to bound shell x cosine temporaries.
            for j=1:128:size(axes,1)
                a=j:min(j+127,size(axes,1));
                for i=1:32:size(dirs,1)
                    d=i:min(i+31,size(dirs,1));
                    u=dirs(d,:)*axes(a,:)';
                    K=obj.response(p,reshape(min(1,abs(u)).^2,1,[]));
                    Z=reshape(K,[],numel(d),numel(a));
                    Y(d,:)=Y(d,:)+reshape(sum(Z.*reshape(mass(a),1,1,[]),3),[],numel(d))';
                end
            end
        end

        function plan = fitPlan(obj,spec)
            plan=tmix.FitPlan(obj,spec);
        end

        function [p,info] = fitNLLS(obj,S,spec,start,opt)
            % Predict(p) optionally evaluates S0/S2 on the measured gradient grid.
            % Ground-truth starts are for synthetic noise propagation only.
            arguments
                obj; S (:,2) double {mustBeReal,mustBeFinite}; spec struct; start struct
                opt.Predict (1,1) function_handle = @(p) obj.invariants(p)
                opt.Weights (:,2) double {mustBeNonnegative,mustBeFinite} = ones(size(S))
            end
            assert(isequal(size(S),[size(obj.protocol,1) 2],size(opt.Weights)), ...
                'Expected one [S0 S2] row and two weights per shell.');
            P=obj.fitPlan(spec); x=P.encode(start);
            info.start=P.decode(x);
            use=obj.protocol(:,1)>0 & opt.Weights>0;
            info.nFree = numel(x);
            info.nObservations = nnz(use);
            assert(isempty(x) || any(use,'all'),'Include a weighted diffusion measurement.');
            residual=@(x) weightedResidual(opt.Predict,P,x,S,opt.Weights,use);
            info.startSSE=sum(residual(x).^2);
            if isempty(x)
                p=info.start; info.exitflag=1; r=residual(x);
            else
                options=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-11, ...
                    'StepTolerance',1e-9,'OptimalityTolerance',1e-9);
                if info.nObservations < info.nFree
                    options.Algorithm = 'levenberg-marquardt';
                end
                [x,~,r,info.exitflag]=lsqnonlin(residual,x,zeros(size(x)),ones(size(x)),options);
                p=P.decode(x);
            end
            info.sse=sum(r.^2);
        end
    end

    methods (Static)
        function p2 = odfP2(axes,mass)
            [axes,mass]=tmix.odf(axes,mass);
            T=axes'*(axes.*mass);
            p2=sqrt(1.5*sum((T-eye(3)/3).^2,'all'));
        end
    end

    methods (Access=private)
        function K = response(obj,p,u2)
            b=obj.protocol(:,1); K=zeros(numel(b),numel(u2));
            if isfield(p,'fStick'), K=K+p.fStick*exp(-b*(p.Dstick*u2)); end
            if isfield(p,'fECS')
                De=p.DePerp+(p.DePar-p.DePerp)*u2;
                if isfield(p,'fPerm') && p.fPerm>0 && p.fECS>0
                    total=p.fPerm+p.fECS;
                    K=K+total*tmix.exchange(obj.protocol,p.Dperm*u2,De, ...
                        p.fPerm/total,p.tau,obj.exchange,obj.steps);
                else
                    K=K+p.fECS*exp(-b*De);
                    if isfield(p,'fPerm'), K=K+p.fPerm*exp(-b*(p.Dperm*u2)); end
                end
            end
            if isfield(p,'fSphere') && p.fSphere>0
                K=K+p.fSphere*tmix.sphere(obj.protocol,p.radius,p.Dsoma);
            end
            if isfield(p,'fDot'), K=K+p.fDot; end
            if isfield(p,'fWater'), K=K+p.fWater*exp(-b*p.Dwater); end
            K(b==0,:)=1;
        end
    end
end

function r=weightedResidual(predict,plan,x,data,weights,use)
prediction=predict(plan.decode(x));
assert(isequal(size(prediction),size(data)),'Predict must return one [S0 S2] row per shell.');
r=(prediction-data).*sqrt(weights); r=r(use);
end
