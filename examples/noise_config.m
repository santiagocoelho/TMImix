function c=noise_config()
% Edit this configuration, then pass the SAME c to either experiment.
[b,D]=ndgrid([0.8 2.5 6],[16 32 56]);
protocol=[b(:) D(:) repmat(5,numel(b),1)];
c.model=TMImix(protocol,["permeableStick","ecs","sphere"],Exchange="finite");
s=c.model.defaults();
s.fPerm=[0.15 0.65]; s.fECS=[0.15 0.75]; s.fSphere=[0.04 0.32];
s.tau=[12 60]; s.radius=[4 10]; s.p2=[0.25 0.85];
c.prior=s;                 % training distribution and free fit bounds
c.testPrior=s;             % independent ground-truth distribution
c.fixed=struct('radius',7); % any active fields; struct() disables overrides
c.mode="matched";          % "matched", "oracle", or "mismatched"
c.snr=[30 60];             % normalized b=0 amplitude / directional noise SD
c.noise="gaussian";        % "rician" for magnitude of complex noisy signals
c.nTissue=8; c.nRep=5; c.nTrain=1400; % nTrain is per SNR
c.seed=712; c.degree=2; c.ridge=1e-4;
c.nDirs=64; c.L=6;         % nuisance l=4,6 only in directional SH projection
end
