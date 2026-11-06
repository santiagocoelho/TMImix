function metrics = exchange_bias_precision(fits)
% Summarize all trials and express biases relative to repeat-measurement SD.
cases = string(fieldnames(fits))';
summaries = cell(numel(cases),2);
for k = 1:numel(cases)
    name = cases(k); result = fits.(name);
    cfg = result.nlls.config;
    fixed = isscalar(cfg.prior.DePar) && isscalar(cfg.prior.DePerp);
    methods = ["nlls","polynomial"];
    for m = 1:2
        t = result.(methods(m)).metrics;
        t.parameter(t.parameter=="fPerm") = "fStick";
        t.parameter(t.parameter=="Dperm") = "Dstick";
        t.fixed = ismember(t.parameter,["DePar","DePerp"]) & fixed;
        t.fixedValue = nan(height(t),1);
        if fixed
            t.fixedValue(t.parameter=="DePar") = cfg.prior.DePar;
            t.fixedValue(t.parameter=="DePerp") = cfg.prior.DePerp;
        end
        t.nTissue = repmat(cfg.nTissue,height(t),1);
        t.nRep = repmat(cfg.nRep,height(t),1);
        t.biasOverRepeatSD = t.bias./t.withinTissueSD;
        t.biasOverRepeatSD(t.fixed | t.withinTissueSD==0) = NaN;
        t = addvars(t,repmat(name,height(t),1),repmat(methods(m),height(t),1), ...
            'Before',1,'NewVariableNames',{'assumption','method'});
        summaries{k,m} = t;
    end
end
metrics = vertcat(summaries{:});
ref = metrics(metrics.assumption=="permeable_free",:);
keys = {'method','parameter'};
[~,refRow] = ismember(metrics(:,keys),ref(:,keys));
metrics.biasShift = metrics.bias-ref.bias(refRow);
matchedSD = ref.withinTissueSD(refRow);
metrics.biasShiftOverMatchedSD = metrics.biasShift./matchedSD;
metrics.biasShiftOverMatchedSD(matchedSD==0) = NaN;
end
