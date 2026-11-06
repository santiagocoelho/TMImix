function plot_bias_precision(metrics,outputDir)
% Compare signed bias with the typical repeat-to-repeat SD in physical units.
params = ["fStick","fECS","fSphere","p2","radius","Dstick","DePar","DePerp","Dsoma","tau"];
labels = ["Stick fraction","ECS fraction","Sphere fraction","fODF p_2","Soma radius (µm)", ...
    "Stick D (µm²/ms)","ECS axial D (µm²/ms)","ECS radial D (µm²/ms)","Soma D (µm²/ms)","Exchange time (ms)"];
cases = ["permeable_free","impermeable_free","permeable_fixed_ECS","impermeable_fixed_ECS"];
caseLabels = ["Permeable · free ECS","Impermeable · free ECS","Permeable · fixed ECS","Impermeable · fixed ECS"];
methods = ["nlls","polynomial"]; methodNames = ["NLLS","Polynomial regression"];
colors = [.06 .32 .53; .72 .30 .07]; symbols = ["o","s"];
fig = figure('Name','Bias and repeat precision','Color','white','Units','pixels','Position',[50 55 1480 610]);
fig.Theme = 'light';
layout = tiledlayout(fig,2,5,'TileSpacing','compact','Padding','compact');
title(layout,sprintf('Bias and repeat precision at b_0 SNR %g',metrics.SNR(1)));
subtitle(layout,'Point = mean fit minus truth   ·   Bar = ±1 within-tissue noise SD   ·   All noise draws included');
for j=1:numel(params)
    ax = nexttile(layout); hold(ax,'on'); title(ax,labels(j),'FontSize',11);
    xline(ax,0,'--','Color',[.44 .44 .44],'HandleVisibility','off');
    for m=1:2
        rows = metrics(metrics.parameter==params(j) & metrics.method==methods(m),:);
        [~,idx] = ismember(rows.assumption,cases);
        y = 5-idx + (1.5-m)*.32;
        errorbar(ax,rows.bias,y,rows.withinTissueSD,'horizontal',symbols(m), ...
            'LineStyle','none','Color',colors(m,:),'MarkerFaceColor',colors(m,:), ...
            'MarkerSize',4,'LineWidth',1.2,'DisplayName',methodNames(m));
    end
    panel = metrics(metrics.parameter==params(j),:);
    bound = max(abs(panel.bias)+panel.withinTissueSD);
    if bound==0, bound=1; end
    xlim(ax,[-1.1 1.1]*bound);
    yticks(ax,1:4); ylim(ax,[.5 4.5]); ax.XGrid='on'; ax.Box='on'; ax.FontSize=9;
    if mod(j-1,5)==0, yticklabels(ax,flip(caseLabels)); else, yticklabels(ax,[]); end
    if j==1
        key = legend(ax,'Orientation','horizontal','Box','off'); key.Layout.Tile='south';
    end
end
xlabel(layout,'Estimation error (native units)');
exportgraphics(fig,fullfile(outputDir,'bias_precision.png'),'Resolution',150,'BackgroundColor','white','Padding',12);
if ~usejava('desktop'), close(fig); end
end
