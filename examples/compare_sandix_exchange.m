function [fits,metrics] = compare_sandix_exchange(c,testTau,fixedTau,fixedECS)
% Fit the same isotropic-ECS measurements with free or prescribed exchange.
% An optional ECS diffusivity adds a fit with fixed De and free exchange.
if nargin<4, fixedECS = []; end
validateattributes(c.snr,{'numeric'},{'scalar','real','positive','nonnan'});
validateattributes(c.nRep,{'numeric'},{'scalar','integer','>=',2});
validateattributes(testTau,{'numeric'},{'row','numel',2,'finite','positive','increasing'});
validateattributes(fixedTau,{'numeric'},{'row','nonempty','finite','positive'});
assert(numel(unique(fixedTau))==numel(fixedTau),'List each fixed exchange time once.');
assert(isequal(c.parts,["permeableStick","ecs","sphere"]) && ...
    c.ecs=="isotropic" && c.mode=="free",'Use isotropic ECS and the free-τ baseline.');
assert(numel(c.prior.tau)==2,'The baseline must estimate exchange time within a range.');
assert(testTau(1)>=c.prior.tau(1) && testTau(2)<=c.prior.tau(2), ...
    'The fit bounds must include every sampled exchange time.');
if ~isempty(fixedECS)
    validateattributes(fixedECS,{'numeric'},{'scalar','real','finite','nonnegative'});
    assert(numel(c.prior.De)==2,'The reference fit must estimate ECS diffusivity.');
end
data = sample_exchange_data(c,testTau);
fits = struct();
caseNames = ["sandix_free",string(matlab.lang.makeValidName(compose('tau_fixed_%.12g_ms',fixedTau)))];
if ~isempty(fixedECS), caseNames(end+1) = "ecs_fixed"; end
assert(numel(unique(caseNames))==numel(caseNames),'Fixed exchange times need distinct case names.');
for k = 1:numel(caseNames)
    cfg = c;
    label = "SANDIX-like · free τ";
    if caseNames(k)=="ecs_fixed"
        cfg.prior.De = fixedECS;
        label = sprintf('Fixed ECS D = %g µm²/ms · free τ',fixedECS);
    elseif k>1
        cfg.prior.tau = fixedTau(k-1);
        label = sprintf('SANDIX-fixed · τ = %g ms',cfg.prior.tau);
        if cfg.prior.tau==1e6, label = "SANDI-like · τ = 10^6 ms"; end
    end
    fprintf('\n%s, b0 SNR %g\n',label,c.snr);
    result = fit_exchange_data(cfg,data,testTau,true);
    result.label = string(label);
    fits.(caseNames(k)) = result;
end
metrics = exchange_bias_precision(fits,"sandix_free");
ref = metrics(metrics.assumption=="sandix_free",:);
[~,ix] = ismember(metrics(:,{'method','parameter'}),ref(:,{'method','parameter'}));
metrics.rmseRatio = metrics.RMSE./ref.RMSE(ix);
metrics.precisionRatio = metrics.withinTissueSD./ref.withinTissueSD(ix);
metrics.rmseRatio(ref.RMSE(ix)==0) = NaN;
metrics.precisionRatio(ref.withinTissueSD(ix)==0) = NaN;
end
