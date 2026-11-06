function c = noise_config(parts,protocol)
% Choose compartments and supply one [b, Delta, Ndirs] row per shell.
% c.delta is the PGSE pulse duration (ms), separate from this table.
% A table with named b, Delta, Ndirs columns is also accepted.
if nargin < 1
    parts = ["permeableStick","ecs","sphere"];
end
if nargin < 2
    [b,Delta] = ndgrid([0.8 2.5 6],[16 32 56]);
    protocol = [b(:), Delta(:), repmat(64,numel(b),1)];
end
c.parts = reshape(string(parts),1,[]);
c.protocol = protocol;
c.delta = 5;
c.exchange = "finite";
% Only the selected compartments determine the parameter names and bounds.
schema = TMImix([1 20 5],c.parts);
c.prior = schema.freeBounds(); % optimizer bounds; scalar fields are fixed
c.testPrior = struct();        % override sampling laws for test tissue
c.trainPrior = struct();       % independent laws for polynomial training
c.fixed = struct();
c.mode = "matched";            % empty c.fixed leaves all active parameters free
c.snr = [30 60];               % true b=0 signal / directional noise SD
c.noise = "gaussian";          % or "rician"
c.nTissue = 8; c.nRep = 5; c.nTrain = 1400;
c.seed = 712; c.degree = 2; c.ridge = 1e-4;
c.L = 6;                       % maximum even SH order allowed per shell
end
