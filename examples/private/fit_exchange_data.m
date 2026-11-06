function result = fit_exchange_data(c,data,exchangeTime,useGrid)
result.forwardTruth = data.truth;
result.exchangeTime = exchangeTime;

% Map the two stick parameters into the fit schema. All other truths carry over.
if ismember("stick",c.parts)
    shared = rmfield(data.truth,{'tau','fPerm','Dperm'});
    [shared.fStick] = deal(data.truth.fPerm);
    [shared.Dstick] = deal(data.truth.Dperm);
    data.truth = orderfields(shared,c.prior);
end
fprintf('Fitting %d trials by NLLS, starting at the corresponding ground truths.\n',height(data.trials));
result.nlls = run_noise(c,"nlls",data,useGrid);
fprintf('Training on %d independent tissues from the fit model per SNR.\n',c.nTrain);
result.polynomial = run_noise(c,"polynomial",data);
end
