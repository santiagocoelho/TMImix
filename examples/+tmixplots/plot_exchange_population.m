function files = plot_exchange_population(population,outputDir)
% Each column uses identical x/y limits across SNRs and both estimators.
groups = {["fStick","fECS","fSphere","p2"], ["Dstick","DePar","DePerp","Dsoma","radius"]};
groupNames = ["fractions_p2","diffusivities_radius"];
titles = {["Stick fraction","ECS fraction","Sphere fraction","fODF p_2"], ...
    ["Stick D (µm²/ms)","ECS axial D (µm²/ms)","ECS radial D (µm²/ms)","Soma D (µm²/ms)","Soma radius (µm)"]};
description = "exchanging tissue fitted with impermeable sticks";
if ismember("permeableStick",population.nlls.config.parts)
    groups{1}(1) = "fPerm";
    groups{2}(1) = "Dperm";
    groups{2}(end+1) = "tau";
    titles{2}(end+1) = "Exchange time (ms)";
    groupNames(2) = "diffusivities_radius_tau";
    description = "matched permeable sticks, ECS and spheres";
end
methods = ["nlls","polynomial"]; methodTitles = ["NLLS","Polynomial regression"];
snrs = unique(population.nlls.trials.snr);
files = strings(0,1);
for g = 1:numel(groups)
    params = groups{g};
    limits = zeros(numel(params),2);
    for col = 1:numel(params)
        name = params(col);
        v = [population.nlls.truth.(name),population.nlls.estimate.(name),population.polynomial.estimate.(name)];
        lo = min(v); hi = max(v);
        pad = 0.04*(hi-lo);
        if pad==0, pad=0.05*max(1,abs(lo)); end
        limits(col,:) = [lo-pad hi+pad];
    end
    for m = 1:numel(methods)
        result = population.(methods(m));
        for first = 1:4:numel(snrs)
            rows = snrs(first:min(first+3,end));
            fig = figure('Name',methodTitles(m)+" · "+groupNames(g),'NumberTitle','off', ...
                'Color','white','Units','pixels','Position',[35 45 270*numel(params) 130+235*numel(rows)]);
            fig.Theme = 'light';
            layout = tiledlayout(fig,numel(rows),numel(params),'TileSpacing','compact','Padding','compact');
            title(layout,sprintf('%s   ·   %s',methodTitles(m),description), ...
                'FontSize',15);
            subtitle(layout,sprintf('Ground-truth τ = %g ms   ·   %d shared tissues   ·   noise realizations per tissue at each b_0 SNR = %d', ...
                population.exchangeTime,result.config.nTissue,result.config.nRep),'FontSize',10);
            xlabel(layout,'Ground truth'); ylabel(layout,'Fitted value');
            for row = 1:numel(rows)
                use = result.trials.snr==rows(row);
                flagged = result.trials.flag(use)<=0;
                for col = 1:numel(params)
                    ax = nexttile(layout); name = params(col);
                    x = [result.truth(use).(name)]'; y = [result.estimate(use).(name)]';
                    handles = drawPanel(ax,x,y,flagged,limits(col,:));
                    if row==1, title(ax,titles{g}(col),'FontSize',11); end
                    if col==1, ylabel(ax,sprintf('SNR %g',rows(row))); end
                    if row<numel(rows), ax.XTickLabel = []; end
                    if row==1 && col==1
                        if ~any(result.trials.flag<=0), handles(end)=[]; end
                        key = legend(ax,handles,'Location','southoutside','Orientation','horizontal', ...
                            'Box','off','AutoUpdate','off');
                        key.Layout.Tile = 'south';
                    end
                end
            end
            stem = methods(m)+"_"+groupNames(g);
            if numel(snrs)>4, stem = stem+sprintf('_page_%02d',ceil(first/4)); end
            files(end+1,1) = fullfile(outputDir,stem+".png"); %#ok<AGROW>
            exportgraphics(fig,files(end),'Resolution',150,'BackgroundColor','white','Padding',12);
            if ~usejava('desktop'), close(fig); end
        end
    end
end
end

function handles = drawPanel(ax,x,y,flagged,lim)
hold(ax,'on');
% Quantile bins give comparable support across a nonuniform tissue prior.
edges = unique(quantile(x,linspace(0,1,min(10,max(1,floor(numel(x)/8)))+1)));
bin = ones(size(x));
if numel(edges)>1, bin = discretize(x,edges); end
[~,~,bin] = unique(bin);
xbar = accumarray(bin,x,[],@mean);
ybar = accumarray(bin,y,[],@mean);
err = y-x;
spread = accumarray(bin,err,[],@std);
spread(accumarray(bin,1)<2) = NaN;
if numel(xbar)>1
    band = fill(ax,[xbar;flipud(xbar)],[ybar-spread;flipud(ybar+spread)], ...
        [.96 .71 .45],'FaceAlpha',.48,'EdgeColor','none');
else
    band = errorbar(ax,xbar,ybar,spread,'Color',[.87 .58 .32],'LineWidth',1.5);
end
band.DisplayName = 'Binned mean ± 1 error SD';
dots = scatter(ax,x,y,9,[.17 .36 .50],'filled','MarkerFaceAlpha',.40,'DisplayName','Individual fits');
identity = plot(ax,lim,lim,'--','Color',[.28 .28 .28],'DisplayName','Identity');
trend = plot(ax,xbar,ybar,'-o','Color',[.56 .22 .04],'MarkerFaceColor',[.56 .22 .04], ...
    'LineWidth',1.3,'MarkerSize',3,'DisplayName','Binned mean');
% A NaN marker preserves the legend handle when this SNR has no flagged exits.
failed = plot(ax,[NaN;x(flagged)],[NaN;y(flagged)],'rx', ...
    'LineWidth',1,'MarkerSize',5,'DisplayName','Flagged NLLS exit');
text(ax,.04,.96,sprintf('Bias %+.3g   SD %.3g',mean(err),std(err)), ...
    'Units','normalized','VerticalAlignment','top','FontSize',8,'BackgroundColor','white','Margin',1);
axis(ax,'square'); xlim(ax,lim); ylim(ax,lim);
ax.XGrid='on'; ax.YGrid='on'; ax.GridAlpha=.12; ax.Box='on'; ax.FontSize=9;
handles = [dots trend band identity failed];
end
