function plot_noise_modes(metrics,modes,units,pdfPath,showFigures)
% One PDF page per parameter. Both methods share an axis, and all modes
% share the y scale for a given parameter and metric.
parameters = unique(metrics.parameter,'stable');
labels = struct('fPerm',"Permeable stick fraction",'fECS',"ECS fraction", ...
    'fSphere',"Sphere fraction",'DePar',"Axial ECS diffusivity", ...
    'DePerp',"Radial ECS diffusivity",'Dperm',"Axial diffusivity of permeable sticks", ...
    'tau',"Exchange time",'radius',"Soma radius",'Dsoma',"Intrasomal diffusivity", ...
    'p2',"fODF anisotropy");
snr = unique(metrics.SNR);
margin = max(1,0.06*(snr(end)-snr(1)));
methods = ["NLLS","Polynomial regression"];
modeTitles = ["Free","Matched","Oracle","Mismatched"];
colors = [0.055 0.32 0.53; 0.76 0.31 0.08];
styles = ["-o","--s"];
sizes = [7,4]; % both markers remain visible when two estimates coincide
visible = "off";
if showFigures, visible = "on"; end
for p = 1:numel(parameters)
    name = parameters(p);
    rows = metrics(metrics.parameter==name,:);
    unit = "";
    if isfield(units,name), unit = " ("+units.(name)+")"; end
    fig = figure('Name',labels.(name),'NumberTitle','off','Color','white', ...
        'Units','pixels','Position',[60 70 1120 610],'Visible',visible);
    fig.Theme = "light";
    layout = tiledlayout(fig,2,numel(modes),'TileSpacing','compact','Padding','compact');
    title(layout,labels.(name)+unit,'Interpreter','none','FontSize',17);
    fits = rows(rows.method=="NLLS",:);
    nTotal = sum(fits.n);
    nConverged = nTotal-sum(fits.nonconverged);
    caption = sprintf('Permeable sticks, ECS and spheres   ·   %d trials per SNR and mode   ·   %d/%d NLLS fits converged', ...
        fits.n(1),nConverged,nTotal);
    subtitle(layout,caption,'FontSize',10);
    for row = 1:2
        if row==1
            metric = 'bias'; axisTitle = "Bias";
            bound = finiteScale(rows.bias);
            limits = [-1.1 1.1]*bound;
        else
            metric = 'RMSE'; axisTitle = "RMSE";
            bound = finiteScale(rows.RMSE);
            limits = [-0.025 1.1]*bound; % keep markers at zero fully visible
        end
        for col = 1:numel(modes)
            ax = nexttile(layout);
            hold(ax,'on');
            ax.YGrid='on'; ax.Box='off'; ax.FontSize=10;
            if row==1
                yline(ax,0,':','Color',[.45 .45 .45],'HandleVisibility','off');
            end
            group = rows(rows.mode==modes(col),:);
            unplottable = strings(0,1);
            for j = 1:numel(methods)
                curve = sortrows(group(group.method==methods(j),:),'SNR');
                values = curve.(metric);
                if any(~isfinite(values))
                    % For example, tau=Inf in a mismatched fit gives infinite error.
                    unplottable(end+1,1) = methods(j);
                    values(~isfinite(values)) = NaN;
                end
                plot(ax,curve.SNR,values,styles(j),'Color',colors(j,:), ...
                    'LineWidth',1.6,'MarkerSize',sizes(j),'DisplayName',methods(j));
            end
            if ~isempty(unplottable)
                text(ax,0.04,0.94,["Nonfinite "+axisTitle; unplottable], ...
                    'Units','normalized','VerticalAlignment','top','FontSize',9);
            end
            xlim(ax,[snr(1)-margin,snr(end)+margin]);
            if numel(snr)<=10, xticks(ax,snr); end
            ylim(ax,limits);
            if row==1
                title(ax,modeTitles(col));
                subtitle(ax,split(group.constraint(1),", "),'Interpreter','none','FontSize',9);
                ax.XTickLabel=[];
            end
            if col==1, ylabel(ax,axisTitle); end
            if row==2 && col==1
                key = legend(ax,'Orientation','horizontal','Box','off');
                key.Layout.Tile = 'south';
            end
        end
    end
    xlabel(layout,'b₀ SNR','FontSize',11);
    exportgraphics(fig,pdfPath,'ContentType','vector','BackgroundColor','white','Padding',14,'Append',p>1);
    if ~showFigures, close(fig); end
end
end


function bound = finiteScale(values)
values = abs(values(isfinite(values)));
bound = max(values);
if isempty(bound) || bound==0, bound=1; end
end
