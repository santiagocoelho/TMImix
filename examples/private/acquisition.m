function [m,acq] = acquisition(c)
% Prepare shell-specific gradient grids and SH projectors once per experiment.
if isfield(c,'model')
    % Previously saved configurations contain the physical model and nDirs.
    m = c.model;
    timing = m.protocol;
    counts = c.nDirs;
    if isscalar(counts), counts = repmat(counts,size(timing,1),1); end
    rows = [timing(:,1:2),counts(:)];
    pulse = timing(:,3);
else
    rows = c.protocol;
    if istable(rows)
        assert(all(ismember(["b","Delta","Ndirs"],string(rows.Properties.VariableNames))), ...
            'Protocol table must contain b, Delta, and Ndirs columns.');
        rows = rows{:,{'b','Delta','Ndirs'}};
    end
    validateattributes(rows,{'numeric'},{'real','finite','2d','nonempty','ncols',3});
    pulse = c.delta(:);
    if isscalar(pulse), pulse = repmat(pulse,size(rows,1),1); end
    validateattributes(pulse,{'numeric'},{'real','finite','nonnegative','numel',size(rows,1)});
    assert(all(rows(:,1:2)>=0,'all'),'b and Delta must be nonnegative.');
    timing = [rows(:,1:2),pulse];
    % Preserve supplied PGSE timing at b=0. Older protocols can use an
    % unspecified Delta or a placeholder shorter than the pulse.
    unspecified = rows(:,1)==0 & (rows(:,2)==0 | pulse>rows(:,2));
    timing(unspecified,:) = repmat([0 1 0],nnz(unspecified),1);
    pulse(unspecified) = 0;
    ecs = "axial";
    if isfield(c,'ecs'), ecs = c.ecs; end
    m = TMImix(timing,c.parts,Exchange=c.exchange,ECS=ecs);
end
counts = rows(:,3);
validateattributes(counts,{'numeric'},{'finite','integer','positive'});
validateattributes(c.L,{'numeric'},{'scalar','integer','>=',2});
assert(mod(c.L,2)==0,'The maximum SH order must be even.');
baseline = rows(:,1)==0;
assert(any(~baseline),'Add at least one diffusion-weighted shell.');
assert(all(counts(~baseline)>=6),'Each nonzero-b shell needs at least 6 directions for S2.');
start = cumsum(counts)-counts;
acq.nMeasurements = sum(counts);
acq.baselineRows = find(baseline);
acq.baselineIndices = find(repelem(baseline,counts));
acq.normalizeB0 = ~isempty(acq.baselineIndices) && ~isfield(c,'model');
acq.groups = struct('shells',{},'index',{},'model',{},'dirs',{},'W',{});
orders = zeros(size(counts));
for n = unique(counts(~baseline))'
    shells = find(counts==n & ~baseline);
    ordersForGrid = 2:2:c.L;
    L = max(ordersForGrid((ordersForGrid+1).*(ordersForGrid+2)/2 <= n));
    orders(shells) = L;
    j = (0:n-1)';
    z = (j+0.5)/n;
    phi = pi*(3-sqrt(5))*j;
    dirs = [sqrt(1-z.^2).*cos(phi),sqrt(1-z.^2).*sin(phi),z];
    g.shells = shells;
    g.index = (1:n)' + start(shells)';
    g.dirs = dirs;
    g.W = tmix.projector(dirs,L);
    if numel(shells)==size(rows,1)
        g.model = m;
    else
        g.model = TMImix(timing(shells,:),m.parts,Exchange=m.exchange,ECS=m.ecs, ...
            Steps=m.steps,Order=numel(m.mu));
    end
    acq.groups(end+1) = g;
end
acq.protocol = array2table([rows,pulse,orders],'VariableNames',{'b','Delta','Ndirs','delta','L'});
acq.shellOfMeasurement = repelem((1:numel(counts))',counts);
end
