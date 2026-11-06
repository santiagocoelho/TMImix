function Y = noise_directions(p,acq)
% A single tissue orientation applies to every shell, including different grids.
v = randn(1,3); v = v/norm(v);
p2 = 0; if isfield(p,'p2'), p2 = p.p2; end
Y = ones(acq.nMeasurements,1);
for g = acq.groups
    S = g.model.invariants(p);
    directional = g.model.kernel(p,g.dirs*v');
    Y(g.index) = (1-p2)*S(:,1)' + p2*directional';
end
end
