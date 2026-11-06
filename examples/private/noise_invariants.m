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
S = zeros(height(acq.protocol),2);
for g = acq.groups
    block = reshape(A(g.index),size(g.index));
    S(g.shells,:) = tmix.measure(block,g.W);
end
if ~isempty(acq.baselineRows)
    means = accumarray(acq.shellOfMeasurement(acq.baselineIndices),A(acq.baselineIndices), ...
        [height(acq.protocol),1],@mean);
    S(acq.baselineRows,1) = means(acq.baselineRows);
end
end
