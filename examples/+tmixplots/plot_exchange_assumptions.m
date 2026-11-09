function overviewFile = plot_exchange_assumptions(fits,outputDir,plotRep)
% One overview of the paired fits, plus a PNG per estimator and assumption.
% Pass [] for plotRep to display every noise draw. The default shows draw 1.
if nargin<3, plotRep = 1; end
cases = string(fieldnames(fits))';
tiles = cell(numel(cases),2); % columns are NLLS and polynomial regression
for k = 1:numel(cases)
    subdir = fullfile(outputDir,cases(k));
    if ~isfolder(subdir), mkdir(subdir); end
    files = tmixplots.plot_exchange_population(fits.(cases(k)),subdir,plotRep,true);
    for m = 1:numel(files)
        tiles{k,m} = imread(files(m));
    end
end

% Copy the exported RGB panels at native resolution so zooming preserves detail.
height = max(cellfun(@(im) size(im,1),tiles),[],'all');
width = max(cellfun(@(im) size(im,2),tiles),[],'all');
[nRows,nCols] = size(tiles);
overview = repmat(uint8(255),height*nRows,width*nCols,3);
for k = 1:nRows
    for m = 1:nCols
        tile = tiles{k,m};
        y = (k-1)*height + (1:size(tile,1));
        x = (m-1)*width + (1:size(tile,2));
        overview(y,x,:) = tile;
    end
end
overviewFile = fullfile(outputDir,'all_parameters.png');
imwrite(overview,overviewFile);
end
