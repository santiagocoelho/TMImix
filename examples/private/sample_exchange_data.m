function data = sample_exchange_data(c,exchangeTime)
forward = c;
forward.parts = ["permeableStick","ecs","sphere"];
[model,acq] = acquisition(forward);
laws = c.testPrior;
laws.tau = exchangeTime; % scalar ground truth or uniform interval in ms
draw = tmix.tissueSampler(model,model.freeBounds(),laws);
data = sample_noise_trials(c,acq,model.defaults(),draw);
end
