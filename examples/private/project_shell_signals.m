function S = project_shell_signals(A,acq)
% Apply the acquisition's SH projection to normalized directional signals.
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
