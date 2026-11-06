%% Locate TMImix when running the saved script or an Editor temporary copy
roots = [string(fileparts(fileparts(mfilename('fullpath')))), string(fileparts(which('TMImix'))), ...
    string(pwd), string(fileparts(pwd)), fullfile(string(getenv('HOME')),'Git','TMImix')];
repoRoot = roots(find(roots~="" & isfile(fullfile(roots,'TMImix.m')),1));
assert(~isempty(repoRoot),'Open the saved example in your TMImix checkout or run it by its full path.');
addpath(repoRoot);

%% Generate exact, b0-normalized S0 and S2 for one exchanging tissue
bValues = [0.8 2.5 6];                 % ms/µm²
diffusionTimes = [16 32 56];           % Δ in ms
[b,Delta] = ndgrid(bValues,diffusionTimes);
protocol = [b(:), Delta(:), 5*ones(numel(b),1)];  % [b, Δ, δ], δ = 5 ms
forwardModel = TMImix(protocol,["permeableStick","ecs","sphere"]);
truth = forwardModel.defaults();
truth.tau = 10;                        % ms; set to Inf for a no-exchange control
truth.fPerm = 0.45; truth.fECS = 0.35; truth.fSphere = 0.20;
truth.Dperm = 2; truth.DePar = 1.5; truth.DePerp = 0.7;   % µm²/ms
truth.radius = 7; truth.Dsoma = 2; truth.p2 = 0.6;        % radius in µm
observed = forwardModel.invariants(truth);

%% Fit those same measurements with impermeable sticks, ECS and spheres
fitModel = TMImix(protocol,["stick","ecs","sphere"]);
bounds = fitModel.freeBounds();        % set any field to a scalar to fix it
sharedTruth = rmfield(truth,{'tau','fPerm','Dperm'});
sharedTruth.fStick = truth.fPerm; sharedTruth.Dstick = truth.Dperm;
[nlls,info] = fitModel.fitNLLS(observed,bounds,sharedTruth);

% Train the polynomial inverse exclusively on noiseless, impermeable tissues.
rng(712,'twister');
nTrain = 2400;
trainPrior = struct('fractionAlpha',[1 1 1]);  % uniform over the fraction simplex
draw = tmix.tissueSampler(fitModel,bounds,trainPrior);
trainTruth = repmat(fitModel.defaults(),nTrain,1);
trainSignal = zeros(size(protocol,1),2,nTrain);
for i = 1:nTrain
    trainTruth(i) = draw();
    trainSignal(:,:,i) = fitModel.invariants(trainTruth(i));
end
regressor = tmix.PolyRegressor(fitModel,bounds,trainSignal,trainTruth,Degree=2);
poly = regressor.predict(observed);

%% Compare physical parameters and forward predictions from both fits
names = string(fieldnames(fitModel.defaults()));
values = @(p) arrayfun(@(name) p.(name),names);
parameterFits = table(names,values(sharedTruth),values(nlls),values(poly), ...
    'VariableNames',{'Parameter','GroundTruth','NLLS','Polynomial'});
nllsSignal = fitModel.invariants(nlls);
polySignal = fitModel.invariants(poly);
signalRMSE = array2table(sqrt([mean((nllsSignal-observed).^2,1); mean((polySignal-observed).^2,1)]), ...
    'VariableNames',{'S0','S2'},'RowNames',{'NLLS','Polynomial'});
fprintf('Ground-truth exchange time = %g ms. Fits use impermeable sticks. NLLS exit flag = %d.\n', ...
    truth.tau,info.exitflag);
disp(parameterFits);
disp('Signal RMSE across shells (normalized to b0)'); disp(signalRMSE);

fig = figure('Name','Exchange mismatch','Color','white','Position',[80 80 880 540]);
fig.Theme = 'light';
layout = tiledlayout(fig,2,numel(bValues),'TileSpacing','compact','Padding','compact');
title(layout,sprintf('Exchange time %g ms · impermeable-stick fits',truth.tau));
for k = 1:2
    for j = 1:numel(bValues)
        ax = nexttile(layout);
        ix = protocol(:,1)==bValues(j);
        [times,order] = sort(protocol(ix,2));
        curves = [observed(ix,k),nllsSignal(ix,k),polySignal(ix,k)];
        h = plot(ax,times,curves(order,1),'-ok',times,curves(order,2),'--s', ...
            times,curves(order,3),':^','LineWidth',1.3);
        h(2).Color = [0.055 0.32 0.53];
        h(3).Color = [0.76 0.31 0.08];
        axis(ax,'padded');
        if numel(times)<=6, xticks(ax,unique(times)); end
        if k==1, title(ax,sprintf('b = %g ms/µm²',bValues(j))); end
        if j==1, ylabel(ax,sprintf('S_%d',2*(k-1))); end
        if k==1 && j==1
            key = legend(h,{'Forward truth','NLLS','Polynomial'},'Orientation','horizontal');
            key.Layout.Tile = 'south';
        end
        grid(ax,'on');
    end
end
xlabel(layout,'Diffusion time Δ (ms)');
outputDir = fullfile(repoRoot,'results','exchange_mismatch');
if ~isfolder(outputDir), mkdir(outputDir); end
exportgraphics(fig,fullfile(outputDir,'signal_fits.png'),'Resolution',150);
