function out = run_noise(c,method)
validateattributes(c.nTissue,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.nRep,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.snr,{'numeric'},{'vector','positive','real','nonempty','nonnan'});
assert(numel(unique(c.snr))==numel(c.snr),'List each b=0 SNR once.');
assert(ismember(c.mode,["free","matched","oracle","mismatched"]),'Unknown constraint mode.');
assert(ismember(c.noise,["gaussian","rician"]),'Unknown directional noise model.');
[m,acq] = acquisition(c);
fixed = c.fixed;
if c.mode=="free"
    fixed = struct();
elseif c.mode=="oracle"
    % Only field names matter here. Use a feasible nominal fit specification
    % for initializing the regressor before conditioning it on each tissue.
    if isstruct(fixed)
        keys = string(fieldnames(fixed))';
    else
        keys = reshape(string(fixed),1,[]);
    end
    assert(all(ismember(keys,string(fieldnames(c.prior)))),'Unknown oracle parameter.');
    fullPlan = m.fitPlan(c.prior);
    nominal = fullPlan.decode(0.5*ones(1,fullPlan.dimension));
    fixed = override(struct(),nominal,keys);
end
assert(isstruct(fixed) && isscalar(fixed),'Use a struct of parameter values for this constraint mode.');
keys = string(fieldnames(fixed))';
assert(all(ismember(keys,string(fieldnames(c.prior)))),'Unknown fixed parameter.');
for k = keys, assert(isscalar(fixed.(k)),'Fixed values must be scalar.'); end
fit = override(c.prior,fixed,keys);
fitPlan = m.fitPlan(fit);
testFixed = struct();
if c.mode=="matched", testFixed = fixed; end
testDraw = tmix.tissueSampler(m,c.prior,c.testPrior,testFixed);

old = rng; finish = onCleanup(@()rng(old)); %#ok<NASGU>
rng(c.seed,'twister');
n = c.nTissue*c.nRep*numel(c.snr); N = size(m.protocol,1);
truth = repmat(m.defaults(),n,1); S = zeros(N,2,n);
tissue = zeros(n,1); snr = zeros(n,1); rep = zeros(n,1); b0Scale = ones(n,1); j = 0;
% Finish all held-out measurements before training or fitting either model.
for i = 1:c.nTissue
    p = testDraw();
    clean = simulate(p,acq);
    for a = reshape(c.snr,1,[])
        for h = 1:c.nRep
            j = j+1; truth(j) = p;
            [S(:,:,j),b0Scale(j)] = observe(clean,acq,a,c.noise);
            tissue(j) = i; snr(j) = a; rep(j) = h;
        end
    end
end
phat = truth; flag = nan(n,1); sse = flag; regressors = cell(1,numel(c.snr));
start = truth([]); trainingTruth = truth([]);
if method=="nlls"
    for j = 1:n
        imposed = fit;
        if c.mode=="oracle", imposed = override(fit,truth(j),keys); end
        [phat(j),info] = m.fitNLLS(S(:,:,j),imposed,truth(j));
        flag(j) = info.exitflag; sse(j) = info.sse; start(j,1) = info.start;
    end
else
    validateattributes(c.nTrain,{'numeric'},{'scalar','integer','>=',2});
    rng(c.seed+1,'twister');
    trainLaws = struct();
    if isfield(c,'trainPrior'), trainLaws = c.trainPrior; end
    known = strings(1,0); trainFixed = fixed;
    if c.mode=="oracle", known = keys; trainFixed = struct(); end
    trainDraw = tmix.tissueSampler(m,c.prior,trainLaws,trainFixed);
    trainingTruth = repmat(m.defaults(),c.nTrain,1);
    clean = zeros(acq.nMeasurements,c.nTrain);
    for i = 1:c.nTrain
        trainingTruth(i) = trainDraw();
        clean(:,i) = simulate(trainingTruth(i),acq);
    end
    for a = 1:numel(c.snr)
        X = zeros(N,2,c.nTrain);
        for i = 1:c.nTrain, X(:,:,i) = observe(clean(:,i),acq,c.snr(a),c.noise); end
        R = tmix.PolyRegressor(m,fit,X,trainingTruth,Degree=c.degree,Ridge=c.ridge,Known=known);
        ix = snr==c.snr(a);
        if isempty(known)
            phat(ix) = R.predict(S(:,:,ix));
        else
            supplied = rmfield(truth(ix),setdiff(fieldnames(truth),cellstr(known)));
            phat(ix) = R.predict(S(:,:,ix),supplied);
        end
        regressors{a} = R;
    end
end
out = struct('config',c,'method',method,'protocol',acq.protocol,'truth',truth, ...
    'trainingTruth',trainingTruth,'estimate',phat,'start',start,'observed',S, ...
    'trials',table(tissue,rep,snr,b0Scale,flag,sse),'regressors',{regressors});
out.metrics = metrics(out);
show = ismember(out.metrics.parameter,[keys,fitPlan.freeFractions,fitPlan.names]);
disp(out.metrics(show,:));
if any(flag<=0), fprintf('%d/%d NLLS fits returned a nonpositive exit flag.\n',nnz(flag<=0),n); end
end

function s = override(s,p,keys)
for k = keys, s.(k) = p.(k); end
end

function Y = simulate(p,acq)
% A single tissue orientation applies to every shell, including different grids.
v = randn(1,3); v = v/norm(v);
p2 = 0; if isfield(p,'p2'), p2 = p.p2; end
Y = ones(acq.nMeasurements,1);
for g = acq.groups
    S = g.model.invariants(p);
    directional = g.model.kernel(p,g.dirs*v');
    Y(g.index) = (1-p2)*S(:,1)' + p2*directional';
end
end

function [S,scale] = observe(clean,acq,snr,noise)
% Unit true baseline gives independent directional noise with sigma=1/SNR.
A = clean + randn(size(clean))/snr;
if noise=="rician", A = hypot(A,randn(size(A))/snr); end
scale = 1;
if acq.normalizeB0
    scale = mean(A(acq.baselineIndices));
    assert(scale>0,'Nonpositive measured b=0 mean. Use more baseline repeats or a higher SNR.');
    A = A/scale;
end
S = zeros(height(acq.protocol),2);
for g = acq.groups
    block = reshape(A(g.index),size(g.index));
    S(g.shells,:) = tmix.measure(block,g.W);
end
if ~isempty(acq.baselineRows)
    means = accumarray(acq.shellOfMeasurement(acq.baselineIndices),A(acq.baselineIndices), ...
        [height(acq.protocol),1],@mean);
    S(acq.baselineRows,1) = means(acq.baselineRows);
end
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
