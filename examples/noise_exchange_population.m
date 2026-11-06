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
data = sample_exchange_data(c,exchangeTime);
simulationSNRs = c.snr;
if nargin==3
    selected = data.trials.snr==fitSNR;
    data.truth = data.truth(selected);
    data.trials = data.trials(selected,:);
    data.observed = data.observed(:,:,selected);
    c.snr = fitSNR;
end
result = fit_exchange_data(c,data,exchangeTime,matched);
result.simulationSNRs = simulationSNRs;
end
