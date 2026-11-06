%% Locate TMImix when running the saved script or an Editor temporary copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example in your TMImix checkout or run it by its full path.');
addpath(repoRoot,fullfile(repoRoot,'examples'));

%% Reuse the reference tissues and directions at a single b0 SNR
[c,exchangeTime] = exchange_population_config("permeableStick"); % same pools for synthesis and fit
targetSNR = 1e6;
% Replace the first SNR in the simulation schedule. Keeping the draw count
% per tissue reproduces the reference sweep's tissues, directions and noise seeds.
if ~ismember(targetSNR,c.snr), c.snr(1) = targetSNR; end
fitModel = TMImix([1 10 6],c.parts);
c.mode = "free";
c.prior = fitModel.freeBounds();  % all ten fields, including tau, vary; fractions sum to one
assert(fitModel.fitPlan(c.prior).dimension==9,'Expected nine independent fit parameters.');

%% Fit and plot only the high-SNR observations
outputDir = fullfile(repoRoot,'results','exchange_high_snr_matched');
if ~isfolder(outputDir), mkdir(outputDir); end
population = noise_exchange_population(c,exchangeTime,targetSNR);
startValues = table2array(struct2table(population.nlls.start));
trueValues = table2array(struct2table(population.forwardTruth));
assert(max(abs(startValues-trueValues),[],'all')<1e-12, ...
    'An NLLS start differs from its forward truth. Check the fit bounds.');
save(fullfile(outputDir,'population.mat'),'population','-v7.3');
plotFiles = tmixplots.plot_exchange_population(population,outputDir);
disp(plotFiles);
