function out=noise_polynomial(c)
% Prior-trained SMI-style polynomial regression with held-out noise trials.
% addpath examples; c=noise_config; out=noise_polynomial(c);
if nargin==0, c=noise_config(); end
out=run_noise(c,"polynomial");
end
