function out = run_noise(c,method,data,matchProjection)
if nargin<4, matchProjection = false; end
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
if nargin < 3
    testFixed = struct();
    if c.mode=="matched", testFixed = fixed; end
    testDraw = tmix.tissueSampler(m,c.prior,c.testPrior,testFixed);
    data = sample_noise_trials(c,acq,m.defaults(),testDraw);
end
truth = data.truth; S = data.observed; trials = data.trials;
n = numel(truth); N = size(S,1);
phat = truth; flag = nan(n,1); sse = flag; regressors = cell(1,numel(c.snr));
startSSE = nan(n,1);
start = truth([]); trainingTruth = truth([]);
if method=="nlls"
    workers = 0;
    if isfield(c,'nllsWorkers'), workers = c.nllsWorkers; end
    fibreAxes = zeros(n,3);
    if matchProjection, fibreAxes = trials.fibreAxis; end
    start = truth;
    parfor (j = 1:n,workers)
        if workers==0 && n>=200 && mod(j-1,100)==0, fprintf('NLLS trial %d/%d\n',j,n); end
        imposed = fit;
        if c.mode=="oracle", imposed = override(fit,truth(j),keys); end
        if matchProjection
            % Use the simulated fibre axis and exactly the same finite-direction
            % measurement operator. The optimizer still sees only S0 and S2.
            v = fibreAxes(j,:);
            predict = @(p) project_shell_signals(noise_directions(p,acq,v),acq);
            [phat(j),info] = m.fitNLLS(S(:,:,j),imposed,truth(j),Predict=predict);
        else
            [phat(j),info] = m.fitNLLS(S(:,:,j),imposed,truth(j));
        end
        startSSE(j) = info.startSSE;
        flag(j) = info.exitflag; sse(j) = info.sse; start(j) = info.start;
    end
else
    validateattributes(c.nTrain,{'numeric'},{'scalar','integer','>=',2});
    old = rng; finish = onCleanup(@()rng(old)); %#ok<NASGU>
    rng(c.seed+1,'twister');
    trainLaws = struct();
    if isfield(c,'trainPrior'), trainLaws = c.trainPrior; end
    known = strings(1,0);
    if c.mode=="oracle", known = keys; end
    % Every constant fit value also constrains the training labels, except
    % per-tissue oracle covariates. Population priors can still differ.
    trainFixed = struct();
    for k = string(fieldnames(fit))'
        if isscalar(fit.(k)) && ~ismember(k,known), trainFixed.(k) = fit.(k); end
    end
    trainDraw = tmix.tissueSampler(m,c.prior,trainLaws,trainFixed);
    trainingTruth = repmat(m.defaults(),c.nTrain,1);
    clean = zeros(acq.nMeasurements,c.nTrain);
    for i = 1:c.nTrain
        trainingTruth(i) = trainDraw();
        clean(:,i) = noise_directions(trainingTruth(i),acq);
    end
    for a = 1:numel(c.snr)
        X = zeros(N,2,c.nTrain);
        for i = 1:c.nTrain, X(:,:,i) = noise_invariants(clean(:,i),acq,c.snr(a),c.noise); end
        R = tmix.PolyRegressor(m,fit,X,trainingTruth,Degree=c.degree,Ridge=c.ridge,Known=known);
        ix = trials.snr==c.snr(a);
        if isempty(known)
            phat(ix) = R.predict(S(:,:,ix));
        else
            supplied = rmfield(truth(ix),setdiff(fieldnames(truth),cellstr(known)));
            phat(ix) = R.predict(S(:,:,ix),supplied);
        end
        regressors{a} = R;
    end
end
trials.flag = flag; trials.sse = sse;
if method=="nlls", trials.startSSE = startSSE; end
out = struct('config',c,'method',method,'protocol',acq.protocol,'truth',truth, ...
    'trainingTruth',trainingTruth,'estimate',phat,'start',start,'observed',S, ...
    'trials',trials,'regressors',{regressors});
out.metrics = metrics(out);
show = ismember(out.metrics.parameter,[keys,fitPlan.freeFractions,fitPlan.names]);
disp(out.metrics(show,:));
if any(flag<=0), fprintf('%d/%d NLLS fits returned a nonpositive exit flag.\n',nnz(flag<=0),n); end
end

function s = override(s,p,keys)
for k = keys, s.(k) = p.(k); end
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
