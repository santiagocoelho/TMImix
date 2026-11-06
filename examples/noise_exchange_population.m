function result = noise_exchange_population(c,exchangeTime,fitSNR)
% Sample exchanging tissues once and fit the identical trials with both estimators.
% Optional fitSNR selects a single condition after simulating the full cohort.
validateattributes(exchangeTime,{'numeric'},{'scalar','real','positive','nonnan'});
exchangingParts = ["permeableStick","ecs","sphere"];
matched = isequal(c.parts,exchangingParts);
assert((matched || isequal(c.parts,["stick","ecs","sphere"])) && ...
    c.mode=="free" && ~isfield(c,'model'),'Expected a free exchanging or impermeable three-pool fit.');
assert(ismember(c.noise,["gaussian","rician"]),'Choose Gaussian or Rician noise.');
if nargin==3
    validateattributes(fitSNR,{'numeric'},{'scalar','positive','real','nonnan'});
    assert(ismember(fitSNR,c.snr),'Requested fit SNR is absent from the simulation schedule.');
end
forward = c;
forward.parts = exchangingParts;
[model,acq] = acquisition(forward);
draw = tmix.tissueSampler(model,model.freeBounds(),c.testPrior,struct('tau',exchangeTime));
data = sample_noise_trials(c,acq,model.defaults(),draw);
result.simulationSNRs = c.snr;
if nargin==3
    selected = data.trials.snr==fitSNR;
    data.truth = data.truth(selected);
    data.trials = data.trials(selected,:);
    data.observed = data.observed(:,:,selected);
    c.snr = fitSNR;
end
result.forwardTruth = data.truth;
result.exchangeTime = exchangeTime;

% Map the two stick parameters into the fit schema. All other truths carry over.
if ~matched
    shared = rmfield(data.truth,{'tau','fPerm','Dperm'});
    [shared.fStick] = deal(data.truth.fPerm);
    [shared.Dstick] = deal(data.truth.Dperm);
    data.truth = orderfields(shared,c.prior);
end
fprintf('Fitting %d trials by NLLS, starting at the corresponding ground truths.\n',height(data.trials));
result.nlls = run_noise(c,"nlls",data,matched);
fprintf('Training on %d independent tissues from the fit model per SNR.\n',c.nTrain);
result.polynomial = run_noise(c,"polynomial",data);
end
