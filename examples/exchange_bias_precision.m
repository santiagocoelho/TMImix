function metrics = exchange_bias_precision(fits,reference)
% Summarize all trials and express biases relative to repeat-measurement SD.
if nargin<2, reference = "permeable_free"; end
cases = string(fieldnames(fits))';
summaries = cell(numel(cases),2);
for k = 1:numel(cases)
    name = cases(k); result = fits.(name);
    cfg = result.nlls.config;
    methods = ["nlls","polynomial"];
    for m = 1:2
        t = result.(methods(m)).metrics;
        t.fixed = false(height(t),1);
        t.fixedValue = nan(height(t),1);
        for j = 1:height(t)
            value = cfg.prior.(t.parameter(j));
            if isscalar(value), t.fixed(j)=true; t.fixedValue(j)=value; end
        end
        t.parameter(t.parameter=="fPerm") = "fStick";
        t.parameter(t.parameter=="Dperm") = "Dstick";
        t.nTissue = repmat(cfg.nTissue,height(t),1);
        t.nRep = repmat(cfg.nRep,height(t),1);
        t.biasOverRepeatSD = t.bias./t.withinTissueSD;
        t.biasOverRepeatSD(t.fixed | t.withinTissueSD==0) = NaN;
        t = addvars(t,repmat(name,height(t),1),repmat(methods(m),height(t),1), ...
            'Before',1,'NewVariableNames',{'assumption','method'});
        if isfield(result,'label'), t.assumptionLabel = repmat(result.label,height(t),1); end
        summaries{k,m} = t;
    end
end
metrics = vertcat(summaries{:});
ref = metrics(metrics.assumption==reference,:);
keys = {'method','parameter'};
[~,refRow] = ismember(metrics(:,keys),ref(:,keys));
assert(all(refRow>0),'The reference fit must contain every reported method and parameter.');
metrics.biasShift = metrics.bias-ref.bias(refRow);
matchedSD = ref.withinTissueSD(refRow);
metrics.biasShiftOverMatchedSD = metrics.biasShift./matchedSD;
metrics.biasShiftOverMatchedSD(matchedSD==0) = NaN;
end
