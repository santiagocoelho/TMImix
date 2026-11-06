function result = noise_exchange_mismatch(c,exchangeTime)
% Sample exchanging tissues once and fit the identical trials with both estimators.
validateattributes(exchangeTime,{'numeric'},{'scalar','real','positive','nonnan'});
assert(isequal(c.parts,["stick","ecs","sphere"]) && c.mode=="free" && ~isfield(c,'model'), ...
    'Use noise_config(["stick","ecs","sphere"],protocol) and c.mode="free".');
assert(ismember(c.noise,["gaussian","rician"]),'Choose Gaussian or Rician noise.');
forward = c;
forward.parts = ["permeableStick","ecs","sphere"];
[model,acq] = acquisition(forward);
draw = tmix.tissueSampler(model,model.freeBounds(),c.testPrior,struct('tau',exchangeTime));
data = sample_noise_trials(c,acq,model.defaults(),draw);
result.forwardTruth = data.truth;
result.exchangeTime = exchangeTime;

% Map the two stick parameters into the fit schema. All other truths carry over.
shared = rmfield(data.truth,{'tau','fPerm','Dperm'});
[shared.fStick] = deal(data.truth.fPerm);
[shared.Dstick] = deal(data.truth.Dperm);
data.truth = orderfields(shared,c.prior);
fprintf('Fitting %d trials by NLLS, starting at the corresponding ground truths.\n',height(data.trials));
result.nlls = run_noise(c,"nlls",data);
fprintf('Training on %d independent impermeable tissues per SNR.\n',c.nTrain);
result.polynomial = run_noise(c,"polynomial",data);
end
