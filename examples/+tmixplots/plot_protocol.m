function plot_protocol(protocol,delta,outputDir)
% Plot each prescribed LTE (b, Delta) pair, retaining actual axis spacing.
if istable(protocol), protocol = protocol{:,{'b','Delta','Ndirs'}}; end
b = protocol(:,1); Delta = protocol(:,2); counts = protocol(:,3);
fig = figure('Name','LTE protocol','Color','white','Position',[70 80 800 470]);
fig.Theme = 'light';
ax = axes(fig); hold(ax,'on');
baseline = b==0;
for n = unique(counts(baseline))'
    ix = baseline & counts==n;
    scatter(ax,Delta(ix),b(ix),60,'s','filled','MarkerFaceColor',[.75 .32 .08], ...
        'DisplayName',sprintf('%g × b_0 per Δ',n));
end
nd = unique(counts(~baseline));
colors = lines(numel(nd));
for j = 1:numel(nd)
    ix = ~baseline & counts==nd(j);
    scatter(ax,Delta(ix),b(ix),58,colors(j,:),'filled', ...
        'DisplayName',sprintf('%g directions per shell',nd(j)));
end
axis(ax,'padded'); xticks(ax,unique(Delta)); yticks(ax,unique(b));
ax.XGrid='on'; ax.YGrid='on'; ax.Box='off';
xlabel(ax,'Diffusion time Δ (ms)'); ylabel(ax,'b (ms/µm²)');
title(ax,sprintf('LTE acquisition   ·   δ = %g ms',delta));
subtitle(ax,sprintf('%d shells   ·   %d diffusion-weighted measurements   ·   %d b_0 measurements', ...
    numel(b),sum(counts(~baseline)),sum(counts(baseline))));
legend(ax,'Location','southoutside','Orientation','horizontal','Box','off');
exportgraphics(fig,fullfile(outputDir,'protocol.png'),'Resolution',150,'BackgroundColor','white','Padding',12);
if ~usejava('desktop'), close(fig); end
end
