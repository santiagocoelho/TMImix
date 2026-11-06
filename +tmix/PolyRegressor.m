classdef PolyRegressor
    % Supervised inverse on normalized shell S0/S2, conditioned optionally
    % on known per-observation parameters. Fixed outputs are reconstructed.
    properties (SetAccess=private)
        model
        spec
        known
        center
        scale
        coef
        terms
        degree
        targets
    end
    methods
        function obj=PolyRegressor(model,spec,S,truth,opt)
            % S: shells x 2 x training samples; truth: array of structs.
            arguments
                model; spec struct; S double {mustBeReal,mustBeFinite}; truth struct
                opt.Degree (1,1) double {mustBeMember(opt.Degree,[1 2 3])} = 2
                opt.Ridge (1,1) double {mustBePositive,mustBeFinite} = 1e-4
                opt.Known (1,:) string = strings(1,0)
            end
            obj.model=model; obj.spec=spec; obj.degree=opt.Degree;
            obj.known=opt.Known; P=model.fitPlan(spec);
            assert(all(ismember(obj.known,string(fieldnames(spec)))) && ...
                all(arrayfun(@(k) isscalar(spec.(k)),obj.known)), ...
                'Known fields need scalar placeholders in the fit spec.');
            n=numel(truth); assert(n>=2 && size(S,3)==n,'Mismatch in training sample counts.');
            X=obj.features(S);
            fields=string(fieldnames(spec))';
            obj.targets=fields(arrayfun(@(k) ~isscalar(spec.(k)),fields));
            T=zeros(n,numel(obj.targets));
            for i=1:n
                if ~isempty(obj.known)
                    [P,v]=obj.condition(truth(i)); X(i,end-numel(v)+1:end)=v;
                end
                x=P.encode(truth(i)); truthCheck=P.decode(x);
                for k=string(fieldnames(spec))'
                    a=truthCheck.(k); b=truth(i).(k);
                    assert(isequal(a,b) || abs(a-b)<1e-9, ...
                        'Training truth for %s is inconsistent with the fit prior.',k);
                end
                for j=1:numel(obj.targets)
                    k=obj.targets(j); bounds=spec.(k);
                    T(i,j)=(truth(i).(k)-bounds(1))/diff(bounds);
                end
            end
            obj.center=mean(X,1); obj.scale=std(X,0,1);
            obj.scale(obj.scale<1e-12)=1;
            d=size(X,2); terms=cell(1,obj.degree);
            for k=1:obj.degree
                % combinations with repetition, including all cross terms
                terms{k}=nchoosek(1:d+k-1,k)-(0:k-1);
            end
            obj.terms=terms; F=obj.design(X);
            % Ridge on standardized-input polynomial terms, unpenalized intercept.
            A=(F'*F)/n+diag([0 repmat(opt.Ridge,1,size(F,2)-1)]);
            obj.coef=A\((F'*T)/n);
        end

        function p=predict(obj,S,knownValues)
            X=obj.features(S); n=size(X,1); P=obj.model.fitPlan(obj.spec);
            if ~isempty(obj.known)
                assert(nargin>=3 && numel(knownValues)==n,'Supply known values for each observation.');
                for i=1:n
                    for j=1:numel(obj.known)
                        X(i,end-numel(obj.known)+j)=knownValues(i).(obj.known(j));
                    end
                end
            end
            U=obj.design(X)*obj.coef;
            p=repmat(P.decode(zeros(1,P.dimension)),n,1);
            for i=1:n
                if ~isempty(obj.known), P=obj.condition(knownValues(i)); end
                candidate=P.template;
                for j=1:numel(obj.targets)
                    k=obj.targets(j); a=obj.spec.(k);
                    candidate.(k)=a(1)+diff(a)*U(i,j);
                end
                p(i)=P.decode(P.encode(candidate));
            end
        end
    end

    methods (Access=private)
        function [P,v]=condition(obj,p)
            s=obj.spec; v=zeros(1,numel(obj.known));
            for j=1:numel(v)
                k=obj.known(j); v(j)=p.(k); s.(k)=v(j);
            end
            assert(all(isfinite(v)),'Per-observation regression covariates must be finite.');
            P=obj.model.fitPlan(s);
        end
        function X=features(obj,S)
            ns=size(obj.model.protocol,1);
            assert(size(S,1)==ns && size(S,2)==2 && all(isfinite(S),'all'), ...
                'Expected one [S0 S2] row per training-protocol shell.');
            X=reshape(S(obj.model.protocol(:,1)>0,:,:),[],size(S,3))';
            X=[X zeros(size(X,1),numel(obj.known))];
        end
        function F=design(obj,X)
            Z=(X-obj.center)./obj.scale;
            F=ones(size(X,1),1+sum(cellfun(@(t) size(t,1),obj.terms))); j=1;
            for k=1:obj.degree
                t=obj.terms{k}; ix=j+(1:size(t,1)); A=ones(size(X,1),numel(ix));
                for a=1:k, A=A.*Z(:,t(:,a)); end
                F(:,ix)=A; j=ix(end);
            end
        end
    end
end
