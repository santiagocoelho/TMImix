function plot_bias_precision(metrics,outputDir,metric)
% Compare signed bias with the typical repeat-to-repeat SD in physical units.
if nargin<3, metric = "bias"; end
assert(ismember(metric,["bias","RMSE","withinTissueSD"]),'Unknown summary metric.');
params = ["fStick","fECS","fSphere","p2","radius","Dstick","De","DePar","DePerp","Dsoma","tau"];
labels = ["Stick fraction","ECS fraction","Sphere fraction","fODF p_2","Soma radius (µm)", ...
    "Stick D (µm²/ms)","ECS D (µm²/ms)","ECS axial D (µm²/ms)","ECS radial D (µm²/ms)", ...
    "Soma D (µm²/ms)","Exchange time (ms)"];
keep = ismember(params,metrics.parameter); params=params(keep); labels=labels(keep);
cases = unique(metrics.assumption,'stable')';
caseLabels = replace(cases,"_"," ");
if ismember("assumptionLabel",string(metrics.Properties.VariableNames))
    for k=1:numel(cases), caseLabels(k)=metrics.assumptionLabel(find(metrics.assumption==cases(k),1)); end
else
    knownCases = ["permeable_free","impermeable_free","permeable_fixed_ECS","impermeable_fixed_ECS"];
    knownLabels = ["Permeable · free ECS","Impermeable · free ECS","Permeable · fixed ECS","Impermeable · fixed ECS"];
    [found,ix] = ismember(cases,knownCases); caseLabels(found) = knownLabels(ix(found));
end
plotTitle = "Bias and repeat precision"; stem = "bias_precision";
caption = "Point = mean fit minus truth   ·   Bar = ±1 within-tissue noise SD   ·   All noise draws included";
axisLabel = "Estimation error (native units)";
if metric=="RMSE"
    plotTitle = "Root mean squared error"; stem = "rmse";
    caption = "Error across tissues and noise draws   ·   Smaller values indicate greater accuracy";
    axisLabel = "RMSE (native units)";
elseif metric=="withinTissueSD"
    plotTitle = "Repeat precision"; stem = "repeat_precision";
    caption = "Within-tissue noise SD   ·   Fixed outputs have zero repeat variation";
    axisLabel = "Within-tissue SD (native units)";
end
methods = ["nlls","polynomial"]; methodNames = ["NLLS","Polynomial regression"];
colors = [.06 .32 .53; .72 .30 .07]; symbols = ["o","s"];
nCols = 5; if ismember("De",params), nCols = 3; end
nRows = ceil(numel(params)/nCols);
fig = figure('Name',plotTitle,'Color','white','Units','pixels', ...
    'Position',[50 55 330*nCols 160+(130+25*numel(cases))*nRows]);
fig.Theme = 'light';
layout = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
title(layout,sprintf('%s at b_0 SNR %g',plotTitle,metrics.SNR(1)));
subtitle(layout,caption);
for j=1:numel(params)
    ax = nexttile(layout); hold(ax,'on'); title(ax,labels(j),'FontSize',11);
    xline(ax,0,'--','Color',[.44 .44 .44],'HandleVisibility','off');
    for m=1:2
        rows = metrics(metrics.parameter==params(j) & metrics.method==methods(m),:);
        [~,idx] = ismember(rows.assumption,cases);
        y = numel(cases)+1-idx + (1.5-m)*.32;
        spread = zeros(height(rows),1);
        if metric=="bias", spread=rows.withinTissueSD; end
        errorbar(ax,rows.(metric),y,spread,'horizontal',symbols(m), ...
            'LineStyle','none','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
            'MarkerSize',4,'LineWidth',1.2,'DisplayName',methodNames(m));
    end
    panel = metrics(metrics.parameter==params(j),:);
    bound = max(abs(panel.(metric)));
    if metric=="bias", bound=max(abs(panel.bias)+panel.withinTissueSD); end
    if bound==0, bound=1; end
    xlim(ax,[-1.1 1.1]*bound);
    if metric~="bias", xlim(ax,[-.05 1.1]*bound); end
    yticks(ax,1:numel(cases)); ylim(ax,[.5 numel(cases)+.5]); ax.XGrid='on'; ax.Box='on'; ax.FontSize=9;
    if mod(j-1,nCols)==0, yticklabels(ax,flip(caseLabels)); else, yticklabels(ax,[]); end
    if j==1
        key = legend(ax,'Orientation','horizontal','Box','off'); key.Layout.Tile='south';
    end
end
xlabel(layout,axisLabel);
exportgraphics(fig,fullfile(outputDir,stem+".png"),'Resolution',150,'BackgroundColor','white','Padding',12);
if ~usejava('desktop'), close(fig); end
end
