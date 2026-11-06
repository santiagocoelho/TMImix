%% Locate TMImix when running the saved script or an Editor temporary copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example in your TMImix checkout or run it by its full path.');
addpath(repoRoot,fullfile(repoRoot,'examples'));

%% Acquisition and population settings
bValues = [0 1 2 5 10 15];          % ms/µm²
diffusionTimes = 10:10:80;          % Δ in ms
pulseWidth = 6;                     % δ in ms
countsPerB = [1 64 64 64 64 64];   % b0 repeats or gradient directions at each Δ
assert(numel(countsPerB)==numel(bValues),'Set one measurement count per b-value.');
[b,Delta] = ndgrid(bValues,diffusionTimes);
protocol = [b(:), Delta(:), repmat(countsPerB(:),numel(diffusionTimes),1)];
c = noise_config(["stick","ecs","sphere"],protocol);  % compartments in both fits
c.delta = pulseWidth; c.L = 8; c.mode = "free";
exchangeTime = 10;                  % ms, fixed in all forward tissues
c.snr = [50 100 250 500];           % per-measurement b0 SNR
c.nTissue = 300; c.nRep = 1;        % same tissues across SNRs and methods
c.nTrain = 6000; c.degree = 2;      % independent impermeable training population
c.noise = "gaussian";               % or "rician"
c.seed = 712;
% Each population can use its own ranges, constants, or sampling functions.
% The three fractions sum to one, here with Dirichlet concentrations [4 4 2].
c.testPrior = struct('fractionAlpha',[4 4 2],'Dperm',[1 2.8],'DePar',[1 2.2], ...
    'DePerp',[0.2 0.85],'radius',[5 10],'Dsoma',[1 2.5],'p2',[0.3 0.9]);
c.trainPrior = rmfield(c.testPrior,'Dperm');
c.trainPrior.Dstick = c.testPrior.Dperm; % edit independently of the test population
% c.prior contains the fit bounds. Set a field to a scalar to fix it in both fits.
outputDir = fullfile(repoRoot,'results','exchange_mismatch');
if ~isfolder(outputDir), mkdir(outputDir); end

%% Plot the protocol, then fit paired noisy S0 and S2 measurements
tmixplots.plot_protocol(c.protocol,c.delta,outputDir);
population = noise_exchange_mismatch(c,exchangeTime);
save(fullfile(outputDir,'population.mat'),'population','-v7.3');

%% Plot from the current fits, or load population.mat to redraw later
plotFiles = tmixplots.plot_exchange_population(population,outputDir);
disp(plotFiles);
