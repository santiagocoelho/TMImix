function plot_noise_scatter(comparison,metrics,modes,units,pdfPath,showFigures)
% Plot the paired held-out estimates at their measured coordinates.
parameters = unique(metrics.parameter,'stable');
snrs = unique(metrics.SNR);
constraints = strings(size(modes));
for col = 1:numel(modes)
    constraints(col) = metrics.constraint(find(metrics.mode==modes(col),1));
end
visible = "off";
if showFigures, visible = "on"; end
for p = 1:numel(parameters)
    name = parameters(p);
    label = noise_parameter_title(name,units);
    blocks = cell(size(modes));
    for col = 1:numel(modes)
        pair = comparison.(modes(col));
        t = pair.nlls.trials(:,{'snr','flag'});
        t.mode = repmat(modes(col),height(t),1);
        t.truth = [pair.nlls.truth.(name)]';
        t.nlls = [pair.nlls.estimate.(name)]';
        t.polynomial = [pair.polynomial.estimate.(name)]';
        blocks{col} = t;
    end
    trials = vertcat(blocks{:});
    limits = parameterLimits(trials{:,{'truth','nlls','polynomial'}});
    hasFlagged = any(trials.flag<=0);
    for first = 1:4:numel(snrs)
        pageSnrs = snrs(first:min(first+3,end));
        fig = figure('Name',"Ground truth vs fit | "+label,'NumberTitle','off', ...
            'Color','white','Units','pixels','Position',[50 50 1020 175+200*numel(pageSnrs)], ...
            'Visible',visible);
        fig.Theme = "light";
        grid = tiledlayout(fig,numel(pageSnrs),numel(modes),'TileSpacing','compact','Padding','compact');
        title(grid,label,'Interpreter','none','FontSize',17);
        subtitle(grid,'Each symbol is one tissue and one noise realization. Identical pairs overlap.', ...
            'FontSize',10);
        xlabel(grid,'Ground truth'); ylabel(grid,'Fitted value');
        for row = 1:numel(pageSnrs)
            for col = 1:numel(modes)
                panel = trials(trials.snr==pageSnrs(row) & trials.mode==modes(col),:);
                ax = nexttile(grid);
                drawPanel(ax,panel,limits,hasFlagged);
                if row==1
                    title(ax,upper(extractBefore(modes(col),2))+extractAfter(modes(col),1));
                    subtitle(ax,split(constraints(col),", "),'Interpreter','none','FontSize',9);
                end
                if col==1
                    ylabel(ax,sprintf('b₀ SNR %g',pageSnrs(row)),'FontWeight','bold');
                else
                    ax.YTickLabel = [];
                end
                if row<numel(pageSnrs), ax.XTickLabel = []; end
                if row==1 && col==1
                    key = legend(ax,'Orientation','horizontal','Box','off','AutoUpdate','off');
                    key.Layout.Tile = 'south';
                end
            end
        end
        exportgraphics(fig,pdfPath,'ContentType','vector','BackgroundColor','white', ...
            'Padding',14,'Append',p>1 || first>1);
        if ~showFigures, close(fig); end
    end
end
end

function drawPanel(ax,t,limits,hasFlagged)
blue = [0.055 0.32 0.53]; orange = [0.76 0.31 0.08];
hold(ax,'on');
plot(ax,limits,limits,'--','Color',[.45 .45 .45],'DisplayName','1:1');
nlls = isfinite(t.truth) & isfinite(t.nlls);
poly = isfinite(t.truth) & isfinite(t.polynomial);
ok = nlls & t.flag>0;
scatter(ax,t.truth(ok),t.nlls(ok),36,blue,'o','LineWidth',1,'DisplayName','NLLS');
scatter(ax,t.truth(poly),t.polynomial(poly),12,orange,'s','filled', ...
    'DisplayName','Polynomial regression');
if hasFlagged
    flagged = nlls & t.flag<=0;
    scatter(ax,t.truth(flagged),t.nlls(flagged),42,blue,'x','LineWidth',1.2, ...
        'DisplayName','NLLS (exit flag ≤ 0)');
end
% tau can be fixed to Inf to disable exchange. Count pairs outside finite axes.
if any(~nlls) || any(~poly)
    note = sprintf('Nonfinite pairs\nNLLS %d, polynomial %d',nnz(~nlls),nnz(~poly));
    text(ax,.04,.95,note,'Units','normalized','VerticalAlignment','top','FontSize',8);
end
axis(ax,'square'); xlim(ax,limits); ylim(ax,limits);
ax.Box='on'; ax.FontSize=9;
end

function limits = parameterLimits(values)
values = values(isfinite(values));
if isempty(values), limits = [0 1]; return; end
limits = [min(values) max(values)];
pad = 0.06*diff(limits);
if pad==0, pad = 0.06*max(1,abs(limits(1))); end
limits = limits+[-pad pad];
end
