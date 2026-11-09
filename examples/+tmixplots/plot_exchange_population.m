function files = plot_exchange_population(population,outputDir,plotRep,combined)
% Each column uses identical x/y limits across SNRs and both estimators.
if nargin<3, plotRep = []; end
if nargin<4, combined = false; end
spec = population.nlls.config.prior;
isotropic = isfield(spec,'De');
groups = {["fStick","fECS","fSphere","p2"], ["Dstick","DePar","DePerp","Dsoma","radius"]};
groupNames = ["fractions_p2","diffusivities_radius"];
titles = {["Stick fraction","ECS fraction","Sphere fraction","fODF p_2"], ...
    ["Stick D (µm²/ms)","ECS axial D (µm²/ms)","ECS radial D (µm²/ms)","Soma D (µm²/ms)","Soma radius (µm)"]};
description = "exchanging tissue fitted with impermeable sticks";
if isotropic
    groups{2} = ["Dstick","De","Dsoma","radius"];
    titles{2} = ["Stick D (µm²/ms)","ECS D (µm²/ms)","Soma D (µm²/ms)","Soma radius (µm)"];
end
if ismember("permeableStick",population.nlls.config.parts)
    groups{1}(1) = "fPerm";
    groups{2}(1) = "Dperm";
    groups{2}(end+1) = "tau";
    titles{2}(end+1) = "Exchange time (ms)";
    groupNames(2) = "diffusivities_radius_tau";
    description = "matched permeable sticks, ECS and spheres";
    if ~isotropic && isscalar(spec.DePar), description = "permeable sticks, ECS and spheres"; end
end
show = true(height(population.nlls.trials),1);
singleRep = ~isempty(plotRep);
if singleRep
    show = population.nlls.trials.rep==plotRep;
    assert(any(show),'The selected noise repetition is absent.');
end
if ~isotropic && isscalar(spec.DePar) && isscalar(spec.DePerp)
    description = description+sprintf('   ·   ECS fixed at (%g, %g) µm²/ms',spec.DePar,spec.DePerp);
end
methods = ["nlls","polynomial"]; methodTitles = ["NLLS","Polynomial regression"];
snrs = unique(population.nlls.trials.snr);
if combined
    assert(isscalar(snrs),'The combined scatter figure takes one b0 SNR.');
    % Keep the same panel positions for both stick models, with tau last.
    radius = groups{2}=="radius";
    groups = {[groups{1},groups{2}(radius),groups{2}(~radius)]};
    titles = {[titles{1},titles{2}(radius),titles{2}(~radius)]};
    groupNames = "all_parameters";
    if ismember("permeableStick",population.nlls.config.parts)
        stickLabel = "Permeable sticks";
    else
        stickLabel = "Impermeable sticks";
    end
    if isotropic
        ecsLabel = "isotropic ECS";
    elseif isscalar(spec.DePar) && isscalar(spec.DePerp)
        ecsLabel = sprintf('fixed ECS D (axial %g, radial %g µm²/ms)',spec.DePar,spec.DePerp);
    else
        ecsLabel = "free ECS diffusivities";
    end
    description = stickLabel+"   ·   "+ecsLabel+sprintf("   ·   b_0 SNR %g",snrs);
end
if isfield(population,'label')
    description = population.label;
    if isotropic, description = description+"   ·   isotropic ECS"; end
    description = description+sprintf("   ·   b_0 SNR %g",snrs(1));
end
files = strings(0,1);
for g = 1:numel(groups)
    params = groups{g};
    limits = zeros(numel(params),2);
    for col = 1:numel(params)
        name = params(col);
        v = [population.nlls.truth(show).(name),population.nlls.estimate(show).(name),population.polynomial.estimate(show).(name)];
        lo = min(v); hi = max(v);
        pad = 0.04*(hi-lo);
        if pad==0, pad=0.05*max(1,abs(lo)); end
        limits(col,:) = [lo-pad hi+pad];
    end
    for m = 1:numel(methods)
        result = population.(methods(m));
        for first = 1:4:numel(snrs)
            rows = snrs(first:min(first+3,end));
            nRows = numel(rows); nCols = numel(params);
            if combined, nRows = 2; nCols = 5; end
            if combined && isotropic, nRows = 3; nCols = 3; end
            fig = figure('Name',methodTitles(m)+" · "+groupNames(g),'NumberTitle','off', ...
                'Color','white','Units','pixels','Position',[35 45 270*nCols 130+235*nRows]);
            fig.Theme = 'light';
            layout = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
            title(layout,sprintf('%s   ·   %s',methodTitles(m),description), ...
                'FontSize',15);
            if isscalar(population.exchangeTime)
                tauLabel = sprintf('Ground-truth τ = %g ms',population.exchangeTime);
            else
                tauLabel = sprintf('Ground-truth τ uniform from %g to %g ms',population.exchangeTime);
            end
            caption = tauLabel+sprintf("   ·   %d shared tissues",result.config.nTissue);
            if combined && ~isotropic, caption = "Simulated exchanging sticks + ECS + spheres   ·   "+caption; end
            if singleRep
                caption = caption+sprintf("   ·   noise draw %d of %d shown",plotRep,result.config.nRep);
            else
                caption = caption+sprintf("   ·   noise realizations per tissue at each b_0 SNR = %d",result.config.nRep);
            end
            subtitle(layout,caption,'FontSize',10);
            xlabel(layout,'Ground truth'); ylabel(layout,'Fitted value');
            for row = 1:numel(rows)
                use = show & result.trials.snr==rows(row);
                flagged = result.trials.flag(use)<=0;
                for col = 1:numel(params)
                    ax = nexttile(layout); name = params(col);
                    x = [result.truth(use).(name)]'; y = [result.estimate(use).(name)]';
                    xlimits = limits(col,:);
                    if name=="tau" && isscalar(spec.tau) && spec.tau>100*max(x)
                        pad = .04*max(1,range(x));
                        xlimits = [min(x)-pad,max(x)+pad];
                    end
                    handles = drawPanel(ax,x,y,flagged,limits(col,:),xlimits);
                    if row==1, title(ax,titles{g}(col),'FontSize',11); end
                    if col==1 && ~combined, ylabel(ax,sprintf('SNR %g',rows(row))); end
                    if row<numel(rows), ax.XTickLabel = []; end
                    if row==1 && col==1
                        if ~any(result.trials.flag(show)<=0), handles(end)=[]; end
                        key = legend(ax,handles,'Location','southoutside','Orientation','horizontal', ...
                            'Box','off','AutoUpdate','off');
                        key.Layout.Tile = 'south';
                    end
                end
            end
            if combined && ~ismember("tau",params)
                ax = nexttile(layout); axis(ax,'off');
                text(ax,.5,.5,{'Exchange time','Fit assumes zero exchange'},'Units','normalized', ...
                    'HorizontalAlignment','center','Color',[.34 .34 .34],'FontSize',11);
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

function handles = drawPanel(ax,x,y,flagged,lim,xlimits)
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
identity = plot(ax,xlimits,xlimits,'--','Color',[.28 .28 .28],'DisplayName','Identity');
trend = plot(ax,xbar,ybar,'-o','Color',[.56 .22 .04],'MarkerFaceColor',[.56 .22 .04], ...
    'LineWidth',1.3,'MarkerSize',3,'DisplayName','Binned mean');
% A NaN marker preserves the legend handle when this SNR has no flagged exits.
failed = plot(ax,[NaN;x(flagged)],[NaN;y(flagged)],'rx', ...
    'LineWidth',1,'MarkerSize',5,'DisplayName','Flagged NLLS exit');
textY=.96; alignment='top';
if range(y)==0 && y(1)>mean(lim), textY=.08; alignment='bottom'; end
text(ax,.04,textY,sprintf('Bias %+.3g   SD %.3g',mean(err),std(err)), ...
    'Units','normalized','VerticalAlignment',alignment,'FontSize',8,'BackgroundColor','white','Margin',1);
axis(ax,'square'); xlim(ax,xlimits); ylim(ax,lim);
ax.XGrid='on'; ax.YGrid='on'; ax.GridAlpha=.12; ax.Box='on'; ax.FontSize=9;
handles = [dots trend band identity failed];
end
