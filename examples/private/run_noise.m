function out=run_noise(c,method)
validateattributes(c.nTissue,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.nRep,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.snr,{'numeric'},{'vector','positive','real','nonempty','nonnan'});
assert(ismember(c.mode,["matched","oracle","mismatched"]),'Unknown constraint mode.');
assert(ismember(c.noise,["gaussian","rician"]),'Unknown directional noise model.');
assert(isscalar(c.fixed) && all(ismember(fieldnames(c.fixed),fieldnames(c.prior))), ...
    'Unknown fixed parameter.');
keys=string(fieldnames(c.fixed))'; m=c.model;
for k=keys, assert(isscalar(c.fixed.(k)),'Fixed values must be scalar.'); end
fit=override(c.prior,c.fixed,keys);
test=c.testPrior; train=fit;
if c.mode=="matched", test=override(test,c.fixed,keys); end
if c.mode=="oracle", train=c.prior; end
fitPlan=m.fitPlan(fit); testPlan=m.fitPlan(test); trainPlan=m.fitPlan(train);
old=rng; finish=onCleanup(@()rng(old)); %#ok<NASGU>
rng(c.seed,'twister'); dirs=sphereDirections(c.nDirs); W=tmix.projector(dirs,c.L);
n=c.nTissue*c.nRep*numel(c.snr); N=size(m.protocol,1);
truth=repmat(m.defaults(),n,1); S=zeros(N,2,n);
tissue=zeros(n,1); snr=zeros(n,1); rep=zeros(n,1); j=0;
% Generate the whole held-out set before either estimator trains/fits.
for i=1:c.nTissue
    p=testPlan.decode(rand(1,testPlan.dimension)); clean=simulate(m,p,dirs);
    for a=reshape(c.snr,1,[])
        for h=1:c.nRep
            j=j+1; truth(j)=p; S(:,:,j)=observe(clean,W,a,c.noise);
            tissue(j)=i; snr(j)=a; rep(j)=h;
        end
    end
end
phat=truth; flag=nan(n,1); sse=flag; regressors=cell(1,numel(c.snr));
start=truth([]);
if method=="nlls"
    for j=1:n
        imposed=fit;
        if c.mode=="oracle", imposed=override(fit,truth(j),keys); end
        [phat(j),info]=m.fitNLLS(S(:,:,j),imposed,truth(j));
        flag(j)=info.exitflag; sse(j)=info.sse; start(j,1)=info.start;
    end
else
    rng(c.seed+1,'twister'); known=strings(1,0);
    if c.mode=="oracle", known=keys; end
    trainingTruth=repmat(m.defaults(),c.nTrain,1);
    clean=zeros(c.nDirs,N,c.nTrain);
    for i=1:c.nTrain
        trainingTruth(i)=trainPlan.decode(rand(1,trainPlan.dimension));
        clean(:,:,i)=simulate(m,trainingTruth(i),dirs);
    end
    for a=1:numel(c.snr)
        X=zeros(N,2,c.nTrain);
        for i=1:c.nTrain, X(:,:,i)=observe(clean(:,:,i),W,c.snr(a),c.noise); end
        R=tmix.PolyRegressor(m,fit,X,trainingTruth,Degree=c.degree,Ridge=c.ridge,Known=known);
        ix=snr==c.snr(a);
        % Only explicitly known covariates can reach the predictor.
        if isempty(known), phat(ix)=R.predict(S(:,:,ix));
        else
            supplied=rmfield(truth(ix),setdiff(fieldnames(truth),cellstr(known)));
            phat(ix)=R.predict(S(:,:,ix),supplied);
        end
        regressors{a}=R;
    end
end
out=struct('config',c,'method',method,'truth',truth,'estimate',phat,'start',start, ...
    'observed',S,'trials',table(tissue,rep,snr,flag,sse),'regressors',{regressors});
out.metrics=metrics(out);
show=ismember(out.metrics.parameter,[keys,fitPlan.freeFractions,fitPlan.names]);
disp(out.metrics(show,:));
if any(flag<=0), fprintf('%d/%d NLLS fits returned a nonpositive exit flag.\n',nnz(flag<=0),n); end
end

function s=override(s,p,keys)
for k=keys, s.(k)=p.(k); end
end

function d=sphereDirections(n)
j=(0:n-1)'; z=(j+0.5)/n; phi=pi*(3-sqrt(5))*j;
d=[sqrt(1-z.^2).*cos(phi),sqrt(1-z.^2).*sin(phi),z];
end

function Y=simulate(m,p,dirs)
% Isotropic ODF + an antipodal fibre of mass p2. The isotropic integral is
% exactly K0; no noisy angular discretization of the isotropic part.
s=m.invariants(p); p2=0; if isfield(p,'p2'), p2=p.p2; end
v=randn(1,3); v=v/norm(v);
Y=(1-p2)*s(:,1)'+p2*m.signal(p,dirs,v,1);
end

function S=observe(Y,W,snr,noise)
A=Y+randn(size(Y))/snr;
if noise=="rician", A=hypot(A,randn(size(Y))/snr); end
S=tmix.measure(A,W);
end

function T=metrics(out)
t=out.trials; rows=cell(0,8);
for a=reshape(out.config.snr,1,[])
    ix=t.snr==a; ids=t.tissue(ix);
    for k=string(fieldnames(out.truth))'
        g=[out.truth(ix).(k)]'; p=[out.estimate(ix).(k)]';
        e=p-g; e(isinf(g)&g==p)=0; % fixed tau=Inf (no exchange)
        within=NaN;
        if out.config.nRep>1
            v=arrayfun(@(i) var(e(ids==i)),unique(ids)); within=sqrt(mean(v));
        end
        rows(end+1,:)={a,k,numel(e),mean(e),std(e),within,sqrt(mean(e.^2)),sum(t.flag(ix)<=0)}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'SNR','parameter','n','bias','errorSD','withinTissueSD','RMSE','nonconverged'});
end
