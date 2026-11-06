function out=noise_nlls(c)
% Ground-truth-initialized NLLS noise propagation. Run from the repo root:
% addpath examples; c=noise_config; out=noise_nlls(c);
if nargin==0, c=noise_config(); end
out=run_noise(c,"nlls");
end
