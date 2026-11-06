function [fits,metrics] = compare_exchange_assumptions(c,exchangeTime,fixedECS)
% Compare four inverses on a single set of exchanging-tissue measurements.
validateattributes(c.snr,{'numeric'},{'scalar','real','positive','nonnan'});
validateattributes(c.nRep,{'numeric'},{'scalar','integer','>=',2});
validateattributes(exchangeTime,{'numeric'},{'row','nonempty','real','finite','positive'});
assert(isscalar(exchangeTime) || (numel(exchangeTime)==2 && exchangeTime(1)<exchangeTime(2)), ...
    'Use a fixed exchange time or an increasing [lo hi] interval in ms.');
assert(isequal(c.parts,["permeableStick","ecs","sphere"]) && c.mode=="free", ...
    'Use the permeable-stick population configuration and free mode.');
assert(all(structfun(@(v) numel(v)==2,c.prior)), ...
    'The comparison requires an all-free reference fit.');
assert(all(exchangeTime>=c.prior.tau(1) & exchangeTime<=c.prior.tau(2)), ...
    'Set c.prior.tau wide enough to start NLLS at every sampled exchange time.');
assert(isequal(sort(fieldnames(fixedECS)),{'DePar';'DePerp'}),'Specify both ECS means.');
validateattributes(fixedECS.DePar,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(fixedECS.DePerp,{'numeric'},{'scalar','real','finite','nonnegative'});
data = sample_exchange_data(c,exchangeTime);
cases = ["permeable_free","impermeable_free","permeable_fixed_ECS","impermeable_fixed_ECS"];
for k = 1:numel(cases)
    name = cases(k);
    cfg = c;
    if startsWith(name,"impermeable")
        cfg.parts(1) = "stick";
        schema = TMImix([0 1 0],cfg.parts);
        cfg.prior = orderfields(impermeableFields(cfg.prior),schema.defaults());
        cfg.trainPrior = impermeableFields(cfg.trainPrior);
    end
    fixed = endsWith(name,"fixed_ECS");
    if fixed
        cfg.prior.DePar = fixedECS.DePar;
        cfg.prior.DePerp = fixedECS.DePerp;
    end
    fprintf('\n%s, b0 SNR %g\n',name,c.snr);
    % All four fits predict S0/S2 on the simulated finite gradient grid.
    result = fit_exchange_data(cfg,data,exchangeTime,true);
    fits.(name) = result;
end
metrics = exchange_bias_precision(fits);
end

function s = impermeableFields(s)
% Retain the training laws and physical bounds for every shared parameter.
if isfield(s,'Dperm'), s.Dstick = s.Dperm; s = rmfield(s,'Dperm'); end
if isfield(s,'fPerm'), s.fStick = s.fPerm; s = rmfield(s,'fPerm'); end
if isfield(s,'tau'), s = rmfield(s,'tau'); end
end
