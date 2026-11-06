function [S,scale] = noise_invariants(clean,acq,snr,noise)
% Unit true baseline gives independent directional noise with sigma=1/SNR.
A = clean + randn(size(clean))/snr;
if noise=="rician", A = hypot(A,randn(size(A))/snr); end
scale = 1;
if acq.normalizeB0
    scale = mean(A(acq.baselineIndices));
    assert(scale>0,'Nonpositive measured b=0 mean. Use more baseline repeats or a higher SNR.');
    A = A/scale;
end
S = project_shell_signals(A,acq);
end
