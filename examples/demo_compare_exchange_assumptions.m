%% Locate TMImix when running the saved script or an Editor temporary copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example in your TMImix checkout or run it by its full path.');
addpath(repoRoot,fullfile(repoRoot,'examples'));

%% User settings
targetSNR = 100;                      % per-measurement b0 SNR
c = exchange_population_config("permeableStick");
testTau = [5 40];                    % uniform exchange times across tissues, in ms
c.trainPrior.tau = testTau;          % edit the polynomial training law independently
c.snr = targetSNR;
c.nRep = 4;                           % at least 2 independent draws per tissue
% Defaults: 300 test tissues and 6000 independent polynomial training tissues.
% The default ECS laws are uniform. Edit the two values if using other laws.
fixedECS = struct('DePar',mean(c.testPrior.DePar),'DePerp',mean(c.testPrior.DePerp));
c.nllsWorkers = 4;                    % set to zero to run serially
if c.nllsWorkers>0 && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('Processes',c.nllsWorkers); end
    c.nllsWorkers = min(c.nllsWorkers,pool.NumWorkers);
else
    c.nllsWorkers = 0;
end
outputDir = fullfile(repoRoot,'results','exchange_assumptions',sprintf('snr_%g',targetSNR));
if ~isfolder(outputDir), mkdir(outputDir); end

%% Fit identical noisy data under all four assumptions
[fits,metrics] = compare_exchange_assumptions(c,testTau,fixedECS);
save(fullfile(outputDir,'comparison.mat'),'fits','metrics','fixedECS','testTau','-v7.3');
writetable(metrics,fullfile(outputDir,'bias_precision.csv'));

%% Redraw from fits and metrics without running the estimation again
% Show one paired draw for every tissue, and use all repeats for precision.
plotRep = 1;  % set [] to show every noise draw against ground truth
for name = string(fieldnames(fits))'
    subdir = fullfile(outputDir,name);
    if ~isfolder(subdir), mkdir(subdir); end
    tmixplots.plot_exchange_population(fits.(name),subdir,plotRep);
end
tmixplots.plot_bias_precision(metrics,outputDir);
fprintf('Wrote figures, bias_precision.csv and comparison.mat to %s\n',outputDir);
