%% Locate TMImix when running the saved script or an Editor temporary copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example in your TMImix checkout or run it by its full path.');
addpath(repoRoot,fullfile(repoRoot,'examples'));

%% Acquisition and population settings
% Edit exchange_population_config.m to change the common protocol or tissue priors.
[c,exchangeTime] = exchange_population_config();
outputDir = fullfile(repoRoot,'results','exchange_mismatch');
if ~isfolder(outputDir), mkdir(outputDir); end

%% Plot the protocol, then fit paired noisy S0 and S2 measurements
tmixplots.plot_protocol(c.protocol,c.delta,outputDir);
population = noise_exchange_mismatch(c,exchangeTime);
save(fullfile(outputDir,'population.mat'),'population','-v7.3');

%% Plot from the current fits, or load population.mat to redraw later
plotFiles = tmixplots.plot_exchange_population(population,outputDir);
disp(plotFiles);
