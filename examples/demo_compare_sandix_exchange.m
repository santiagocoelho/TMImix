%% Locate TMImix, including when running a temporary Editor copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example or run it by its full path.');
addpath(repoRoot,fullfile(repoRoot,'examples'));

%% User settings
targetSNR = [50 1e6];                 % b0 SNRs, or a single value
testTau = [5 30];                     % uniform exchange times in ms
fixedTau = [10 1e6];                  % add values, e.g. [5 10 15 20 25 30 1e6]
fixedECS = 1;                        % fourth fit: fixed De, free τ; [] omits it
somaD = 2;                           % fixed in test tissues, training and estimation
p2Range = [0.1 0.5];                 % test/training distribution and fit bounds
bValues = [0 1 2 5 10 15];           % ms/µm²
diffusionTimes = 10:10:80;            % Δ, ms
pulseWidth = 6;                      % δ, ms
nDirections = 30;                    % per nonzero-b shell
nB0 = 1;                            % b0 images per Δ
[b,Delta] = ndgrid(bValues,diffusionTimes);
protocol = [b(:),Delta(:),nDirections*ones(numel(b),1)];
protocol(protocol(:,1)==0,3) = nB0;
c = noise_config(["permeableStick","ecs","sphere"],protocol,"isotropic");
c.delta = pulseWidth; c.L = 6; c.mode = "free";
c.nTissue = 300; c.nRep = 4; c.nTrain = 6000;
c.noise = "gaussian";                % or "rician", applied before S0/S2 extraction
c.degree = 2; c.seed = 712;
% Fractions follow a joint Dirichlet law in [stick, ECS, sphere] order.
% Other fields accept uniform [lo hi] intervals, constants, or @() draws.
c.testPrior = struct('fractionAlpha',[4 4 2],'Dperm',[1 2.8],'De',[0.5 1.5], ...
    'radius',[5 10],'Dsoma',somaD,'p2',p2Range);
c.trainPrior = c.testPrior;          % edit polynomial training laws independently
c.trainPrior.tau = testTau;
c.prior.tau = testTau;               % free-τ fit bounds, editable independently
c.prior.Dsoma = somaD;
c.prior.p2 = p2Range;
% tau = 1/(k_stick_to_ECS + k_ECS_to_stick), following TMImix's convention.
c.nllsWorkers = 4;                   % zero selects serial NLLS
if c.nllsWorkers>0 && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('Processes',c.nllsWorkers); end
    c.nllsWorkers = min(c.nllsWorkers,pool.NumWorkers);
else
    c.nllsWorkers = 0;
end

%% Fit and plot each SNR using the same seeded tissues and standardized noise draws
validateattributes(targetSNR,{'numeric'},{'vector','nonempty','real','positive','nonnan'});
plotRep = 1;                         % [] displays every noise draw against truth
for snr = reshape(targetSNR,1,[])
    c.snr = snr;
    outputDir = fullfile(repoRoot,'results','sandix_exchange_fixed_soma',sprintf('snr_%.12g',snr));
    if ~isfolder(outputDir), mkdir(outputDir); end
    [fits,metrics] = compare_sandix_exchange(c,testTau,fixedTau,fixedECS);
    save(fullfile(outputDir,'comparison.mat'),'fits','metrics','c','testTau','fixedTau','fixedECS','-v7.3');
    writetable(metrics,fullfile(outputDir,'bias_precision.csv'));
    tmixplots.plot_exchange_assumptions(fits,outputDir,plotRep);
    tmixplots.plot_bias_precision(metrics,outputDir);
    tmixplots.plot_bias_precision(metrics,outputDir,"RMSE");
    tmixplots.plot_bias_precision(metrics,outputDir,"withinTissueSD");
    tmixplots.plot_protocol(c.protocol,c.delta,outputDir);
    fprintf('Wrote figures, bias_precision.csv and comparison.mat to %s\n',outputDir);
end
