function data = sample_noise_trials(c,acq,template,testDraw)
% Generate the paired observations before fitting or training.
validateattributes(c.nTissue,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.nRep,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.snr,{'numeric'},{'vector','positive','real','nonempty','nonnan'});
assert(numel(unique(c.snr))==numel(c.snr),'List each b=0 SNR once.');
old = rng; finish = onCleanup(@()rng(old)); %#ok<NASGU>
rng(c.seed,'twister');
n = c.nTissue*c.nRep*numel(c.snr); N = height(acq.protocol);
fibreAxis = zeros(n,3);
truth = repmat(template,n,1); S = zeros(N,2,n);
tissue = zeros(n,1); snr = zeros(n,1); rep = zeros(n,1); b0Scale = ones(n,1); j = 0;
% Finish all held-out measurements before training or fitting either model.
for i = 1:c.nTissue
    p = testDraw();
    [clean,v] = noise_directions(p,acq);
    for a = reshape(c.snr,1,[])
        for h = 1:c.nRep
            j = j+1; truth(j) = p; fibreAxis(j,:) = v;
            [S(:,:,j),b0Scale(j)] = noise_invariants(clean,acq,a,c.noise);
            tissue(j) = i; snr(j) = a; rep(j) = h;
        end
    end
end
data = struct('truth',truth,'observed',S,'trials',table(tissue,rep,snr,b0Scale,fibreAxis));
end
