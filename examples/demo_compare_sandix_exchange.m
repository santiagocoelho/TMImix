%% Locate TMImix, including when running a temporary Editor copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example or run it by its full path.');
addpath(repoRoot,fullfile(repoRoot,'examples'));

%% User settings
targetSNR = 50;                       % b0 signal / per-image noise SD
testTau = [5 30];                     % uniform exchange times in ms
fixedTau = [10 1e6];                  % add values, e.g. [5 10 15 20 25 30 1e6]
bValues = [0 1 2 5 10 15];           % ms/µm²
diffusionTimes = 10:10:80;            % Δ, ms
pulseWidth = 6;                      % δ, ms
nDirections = 30;                    % per nonzero-b shell
nB0 = 1;                            % b0 images per Δ
[b,Delta] = ndgrid(bValues,diffusionTimes);
protocol = [b(:),Delta(:),nDirections*ones(numel(b),1)];
protocol(protocol(:,1)==0,3) = nB0;
c = noise_config(["permeableStick","ecs","sphere"],protocol,"isotropic");
c.delta = pulseWidth; c.L = 6; c.snr = targetSNR; c.mode = "free";
c.nTissue = 300; c.nRep = 4; c.nTrain = 6000;
c.noise = "gaussian";                % or "rician", applied before S0/S2 extraction
c.degree = 2; c.seed = 712;
% Fractions follow a joint Dirichlet law in [stick, ECS, sphere] order.
% Other fields accept uniform [lo hi] intervals, constants, or @() draws.
c.testPrior = struct('fractionAlpha',[4 4 2],'Dperm',[1 2.8],'De',[0.5 1.5], ...
    'radius',[5 10],'Dsoma',[1 2.5],'p2',[0.3 0.9]);
c.trainPrior = c.testPrior;          % edit polynomial training laws independently
c.trainPrior.tau = testTau;
c.prior.tau = testTau;               % free-τ fit bounds, editable independently
% tau = 1/(k_stick_to_ECS + k_ECS_to_stick), following TMImix's convention.
c.nllsWorkers = 4;                   % zero selects serial NLLS
if c.nllsWorkers>0 && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('Processes',c.nllsWorkers); end
    c.nllsWorkers = min(c.nllsWorkers,pool.NumWorkers);
else
    c.nllsWorkers = 0;
end
outputDir = fullfile(repoRoot,'results','sandix_exchange',sprintf('snr_%g',targetSNR));
if ~isfolder(outputDir), mkdir(outputDir); end

%% One set of exchanging tissues and noise draws for every estimator and fit
[fits,metrics] = compare_sandix_exchange(c,testTau,fixedTau);
save(fullfile(outputDir,'comparison.mat'),'fits','metrics','c','testTau','fixedTau','-v7.3');
writetable(metrics,fullfile(outputDir,'bias_precision.csv'));

%% Redraw from saved fits without running the estimation again
plotRep = 1;                         % [] displays every noise draw against truth
tmixplots.plot_exchange_assumptions(fits,outputDir,plotRep);
tmixplots.plot_bias_precision(metrics,outputDir);
tmixplots.plot_bias_precision(metrics,outputDir,"RMSE");
tmixplots.plot_bias_precision(metrics,outputDir,"withinTissueSD");
tmixplots.plot_protocol(c.protocol,c.delta,outputDir);
fprintf('Wrote figures, bias_precision.csv and comparison.mat to %s\n',outputDir);
