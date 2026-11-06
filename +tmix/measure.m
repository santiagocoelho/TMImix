function S=measure(Y,W)
% Baseline-normalized directional DWI [directions x shells] to [shells x 2].
% For shell-specific sampling, call separately with each shell's projector.
assert(size(W,1)==6 && size(W,2)==size(Y,1) && all(isfinite(Y),'all'), ...
    'Expected finite directional signals matching the projector.');
C=W*Y; S=[abs(C(1,:))'/sqrt(4*pi),vecnorm(C(2:6,:),2,1)'/sqrt(20*pi)];
end
