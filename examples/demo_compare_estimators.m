%% User settings
% From the repository root: run('examples/demo_compare_estimators.m')
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot,fullfile(repoRoot,'examples'));

% Edit these independently. Each field is a physical parameter name.
matchedFixed    = struct('radius',6);   % test tissue and fit, in µm
mismatchedFixed = struct('radius',7);   % fit only, in µm
oracleKnown     = "radius";            % supply the true value for each tissue
% Multiple fields are allowed, e.g. struct('radius',6,'Dsoma',2).
% In oracle mode, use a list such as ["radius","Dsoma"].

parts = ["permeableStick","ecs","sphere"];
diffusionTimes = [16 32 56];
[b,Delta] = ndgrid([0.8 2.5 6],diffusionTimes);
Ndirs = repmat([20;32;48],numel(diffusionTimes),1);
protocol = [0 0 6; b(:),Delta(:),Ndirs];    % [b (ms/µm²), Delta (ms), Ndirs]
config = noise_config(parts,protocol);
config.delta = 5;                           % PGSE pulse width (ms)
config.snr = [20 40 60 80];                 % SNR of one b=0 measurement
config.nTissue = 8; config.nRep = 3;        % held-out tissues and noise repeats
config.nTrain = 1400; config.seed = 712;
config.testPrior.radius = [6 10];           % uniform in µm whenever radius varies
config.testPrior.fractionAlpha = [4 4 2];   % in the order given by parts
config.trainPrior = config.testPrior;       % edit training laws independently here

outputDir = fullfile(repoRoot,'results','compare_estimators');
showFigures = usejava('desktop');           % also export all figures in batch mode

%% Run both estimators on the same measurements within each mode
modes = ["free","matched","oracle","mismatched"];
fixedByMode = struct('free',struct(),'matched',matchedFixed, ...
    'oracle',oracleKnown,'mismatched',mismatchedFixed);
units = struct('DePar',"µm²/ms",'DePerp',"µm²/ms",'Dperm',"µm²/ms", ...
    'Dsoma',"µm²/ms",'radius',"µm",'tau',"ms");
comparison = struct();
methodKeys = ["nlls","polynomial"];
methodNames = ["NLLS","Polynomial regression"];
metricTables = cell(numel(modes),2);
validateattributes(config.snr,{'numeric'},{'vector','real','finite','positive','nonempty'});
for i = 1:numel(modes)
    mode = modes(i);
    c = config;
    c.mode = mode;
    c.fixed = fixedByMode.(mode);
    fprintf('\nRunning %s mode\n',mode);
    nlls = noise_nlls(c);
    poly = noise_polynomial(c);
    assert(isequal(nlls.truth,poly.truth) && isequal(nlls.observed,poly.observed), ...
        'The two estimators must see identical test tissues and noise.');
    comparison.(mode) = struct('nlls',nlls,'polynomial',poly);
    constraint = constraintText(mode,c.fixed,units);
    for column = 1:numel(methodKeys)
        t = comparison.(mode).(methodKeys(column)).metrics;
        n = height(t);
        metricTables{i,column} = addvars(t,repmat(mode,n,1),repmat(methodNames(column),n,1), ...
            repmat(constraint,n,1),'Before',1,'NewVariableNames',{'mode','method','constraint'});
    end
end
metrics = vertcat(metricTables{:});
if ~isfolder(outputDir), mkdir(outputDir); end
writetable(metrics,fullfile(outputDir,'metrics.csv'));

%% Plot signed bias (estimate minus truth) and RMSE, in native units
% Run this section again to redraw from the metrics table in the workspace.
pdfPath = fullfile(outputDir,'bias_rmse_by_parameter.pdf');
plot_noise_modes(metrics,modes,units,pdfPath,showFigures);
fprintf('Wrote %s and metrics.csv\n',pdfPath);

%% Scatter plots of ground truth versus fitted values
% Rerun this section with comparison in the workspace to reuse the same fits.
scatterPath = fullfile(outputDir,'ground_truth_vs_fit.pdf');
plot_noise_scatter(comparison,metrics,modes,units,scatterPath,showFigures);
fprintf('Wrote %s\n',scatterPath);

function label = constraintText(mode,fixed,units)
if mode=="oracle"
    names = reshape(string(fixed),1,[]);
    entries = names + " = truth";
else
    names = reshape(string(fieldnames(fixed)),1,[]);
    entries = strings(size(names));
    for j = 1:numel(names)
        name = names(j);
        entries(j) = sprintf('%s = %.6g',name,fixed.(name));
        if isfield(units,name), entries(j) = entries(j)+" "+units.(name); end
    end
end
if isempty(entries), label = "All free"; else, label = join(entries,", "); end
end
