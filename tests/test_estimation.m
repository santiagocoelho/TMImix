function tests=test_estimation
tests=functiontests(localfunctions);
end

function setupOnce(t)
root=fileparts(fileparts(mfilename('fullpath')));
t.applyFixture(matlab.unittest.fixtures.PathFixture(root));
t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'examples')));
end

function testFixedAndBoundedSimplex(t)
m=TMImix([1 20 5;2 40 8],["stick","ecs","dot","water"]);
s=m.defaults(); s.fStick=[.1 .4]; s.fECS=.33; s.fDot=[.12 .45]; s.fWater=[.07 .3];
s.Dstick=[1 3]; P=m.fitPlan(s);
for a=[0 .2 .7 1]
    for b=[0 .5 1]
        p=P.decode([a b .7]); F=[p.fStick p.fECS p.fDot p.fWater];
        verifyEqual(t,sum(F),1,'AbsTol',1e-14); verifyEqual(t,p.fECS,.33);
        verifyGreaterThanOrEqual(t,F,[.1 .33 .12 .07]-1e-14);
        verifyLessThanOrEqual(t,F,[.4 .33 .45 .3]+1e-14);
        q=P.decode(P.encode(p)); verifyEqual(t,cell2mat(struct2cell(q)),cell2mat(struct2cell(p)),'AbsTol',1e-14);
    end
end
end

function testNoiselessNLLS(t)
m=TMImix([.5 20 5;1.5 20 5;3 50 5;6 50 5],["stick","ecs"]);
p=m.defaults(); p.fStick=.42; p.fECS=.58; p.p2=.65;
s=p; s.fStick=[.1 .9]; s.fECS=[.1 .9]; s.p2=[.1 .9];
init=p; init.fStick=.65; init.fECS=.35; init.p2=.4;
[h,info]=m.fitNLLS(m.invariants(p),s,init);
verifyGreaterThan(t,info.exitflag,0); verifyEqual(t,[h.fStick h.p2],[p.fStick p.p2],'AbsTol',1e-6);
verifyEqual(t,h.DePar,p.DePar);
end

function testThreeConstraintModes(t)
c=noise_config(); c.nTissue=3; c.nRep=2; c.nTrain=55; c.snr=45;
c.fixed=struct('radius',6,'fPerm',.35);
for mode=["matched","mismatched","oracle"]
    c.mode=mode; a=noise_nlls(c); b=noise_polynomial(c);
    verifyEqual(t,a.truth,b.truth); verifyEqual(t,a.observed,b.observed);
    r=[a.truth.radius]; f=[a.truth.fPerm];
    if mode=="matched"
        verifyEqual(t,r,6+0*r); verifyEqual(t,f,.35+0*f);
    else
        verifyGreaterThan(t,range(r),.1); verifyGreaterThan(t,range(f),.01);
    end
    if mode=="mismatched", signals=a.observed; end
    if mode=="oracle", verifyEqual(t,a.observed,signals); end
    for result={a,b}
        h=result{1}.estimate;
        if mode=="oracle"
            verifyEqual(t,[h.radius],r); verifyEqual(t,[h.fPerm],f);
        else
            verifyEqual(t,[h.radius],6+0*r); verifyEqual(t,[h.fPerm],.35+0*f);
        end
        verifyEqual(t,[h.fSphere]+[h.fECS]+[h.fPerm],ones(size(r)),'AbsTol',1e-12);
    end
    if mode~="mismatched"
        verifyEqual(t,[a.start.fECS],[a.truth.fECS],'AbsTol',1e-12);
        verifyEqual(t,[a.start.tau],[a.truth.tau],'AbsTol',1e-12);
    end
end
end

function testAllFixedAndIsotropic(t)
m=TMImix([0 20 0;1 30 0],"dot"); p=m.defaults();
S=m.invariants(p); P=m.fitPlan(p);
verifyEqual(t,P.dimension,0); verifyEqual(t,S,[1 0;1 0],'AbsTol',1e-14); verifyEqual(t,m.fitNLLS(S,p,p),p);
reg=tmix.PolyRegressor(m,p,cat(3,S,S),[p p],Degree=3);
verifyEqual(t,reg.predict(S),p);
end

function testPolynomialPredictsUnseenP2(t)
m=TMImix([.7 20 4;2 40 4],"stick"); p=m.defaults(); s=p; s.p2=[.1 .9];
pt=repmat(p,40,1); S=zeros(2,2,40);
for i=1:numel(pt)
    pt(i).p2=.1+.8*(i-1)/39;
    S(:,:,i)=m.invariants(pt(i));
end
r=tmix.PolyRegressor(m,s,S,pt,Degree=2,Ridge=1e-8);
p.p2=.526; h=r.predict(m.invariants(p));
verifyEqual(t,h.p2,p.p2,'AbsTol',1e-6); verifyEqual(t,h.Dstick,p.Dstick);
end

function testAllSixFreeAndMixedProtocol(t)
parts = ["stick","ecs","permeableStick","sphere","dot","water"];
c = noise_config(parts);
n = size(c.protocol,1);
c.protocol(:,3) = repmat([6;15;32],n/3,1);
c.protocol = array2table([0 0 3; c.protocol; 0 0 5], ...
    'VariableNames',{'b','Delta','Ndirs'});
c.delta = [0;repmat([3;4;5],n/3,1);0];
c.nTissue=1; c.nRep=1; c.nTrain=22; c.snr=45; c.mode="free";
c.testPrior.fractionAlpha=2*ones(1,6);
a=noise_nlls(c); b=noise_polynomial(c);
verifyEqual(t,a.truth,b.truth); verifyEqual(t,a.observed,b.observed);
verifyEqual(t,a.protocol.Ndirs,c.protocol.Ndirs);
verifyEqual(t,a.protocol.L,[0; repmat([2;4;6],n/3,1); 0]);
verifyEqual(t,a.protocol.delta,c.delta);
verifyEqual(t,sum(a.observed([1 end],1).*[3;5])/8,1,'AbsTol',1e-14);
verifyEqual(t,a.observed([1 end],2),[0;0]);
m=TMImix([1 20 5],parts);
verifyEqual(t,m.fitPlan(c.prior).dimension,14); % 15 physical fields, one sum constraint
verifyTrue(t,all(structfun(@(v) numel(v)==2,c.prior)));
verifyEqual(t,cell2mat(struct2cell(a.start)),cell2mat(struct2cell(a.truth)),'AbsTol',1e-12);
verifyTrue(t,all(isfinite(cell2mat(struct2cell(b.estimate)))));
end

function testDistinctPopulationsAndConstraints(t)
c = noise_config(["sphere","dot"],[0 0 3;1 15 6;2 30 15;4 45 32]);
c.nTissue=7; c.nRep=1; c.nTrain=40; c.snr=40;
c.testPrior.radius = [6 10];
c.trainPrior.radius = @() 4+rand^2;
c.testPrior.fDot = @() .1+.1*rand;
c.trainPrior.fDot = [.32 .41];
for mode = ["matched","mismatched","oracle"]
    c.mode=mode;
    switch mode
        case "matched", c.fixed=struct('radius',6);
        case "mismatched", c.fixed=struct('radius',7);
        case "oracle", c.fixed="radius";
    end
    a=noise_nlls(c); b=noise_polynomial(c);
    verifyEqual(t,a.observed,b.observed);
    rt=[a.truth.radius]; rr=[b.trainingTruth.radius];
    if mode=="matched"
        verifyEqual(t,rt,6*ones(size(rt))); verifyEqual(t,rr,6*ones(size(rr)));
        expected=6*ones(size(rt));
    else
        verifyTrue(t,all(rt>=6 & rt<=10)); verifyGreaterThan(t,range(rt),.5);
        if mode=="oracle"
            verifyEqual(t,a.observed,saved);
            verifyTrue(t,all(rr>=4 & rr<=5)); expected=rt;
        else
            saved=a.observed;
            verifyEqual(t,rr,7*ones(size(rr))); expected=7*ones(size(rt));
        end
    end
    verifyEqual(t,[a.estimate.radius],expected); verifyEqual(t,[b.estimate.radius],expected);
    verifyTrue(t,all([b.trainingTruth.fDot]>=.32 & [b.trainingTruth.fDot]<=.41));
    verifyTrue(t,all([a.truth.fDot]>=.1 & [a.truth.fDot]<=.2));
end
end

function testB0NoiseAndDirectionCount(t)
% With a pure dot, L=0 coefficient noise and baseline noise are known.
c=noise_config("dot",[0 0 5;1 20 15;1 40 60]);
c.nTissue=1; c.nRep=1500; c.nTrain=12; c.snr=20;
r=noise_polynomial(c); sd=1/c.snr;
verifyEqual(t,std(r.trials.b0Scale),sd/sqrt(5),'RelTol',.12);
% Undo the shared normalization before comparing the two shell means.
raw=squeeze(r.observed(2:3,1,:)).*r.trials.b0Scale';
verifyEqual(t,std(raw(1,:)),sd/sqrt(15),'RelTol',.18);
verifyEqual(t,std(raw(2,:)),sd/sqrt(60),'RelTol',.12);
end

function testUnderDeterminedFit(t)
m = TMImix([1.5 30 4],["stick","ecs"]);
p = m.defaults(); s = m.freeBounds(); obs = m.invariants(p)+[.004 .001];
[q,info] = m.fitNLLS(obs,s,p);
verifyGreaterThan(t,info.nFree,info.nObservations);
verifyGreaterThan(t,info.exitflag,0);
verifyLessThan(t,info.sse,sum((obs-m.invariants(p)).^2));
verifyEqual(t,q.fStick+q.fECS,1,'AbsTol',1e-12);
end

function testPairedExchangeMismatch(t)
protocol = [0 10 1; 1 10 24; 4 10 32; 0 40 2; 1 40 24; 4 40 32];
c = noise_config(["stick","ecs","sphere"],protocol);
c.delta=6; c.mode="free"; c.snr=[Inf 50]; c.nTissue=3; c.nRep=2; c.nTrain=32;
c.testPrior.radius = [5 10]; c.testPrior.fractionAlpha = [4 4 2];
c.prior.radius = 7; c.trainPrior.radius = [5 10];
p = noise_exchange_mismatch(c,10);
a = p.nlls; b = p.polynomial;
verifyEqual(t,a.observed,b.observed);
verifyEqual(t,a.truth,b.truth);
verifyEqual(t,a.start(1).radius,7);
verifyEqual(t,[a.estimate.radius],7*ones(1,12));
verifyEqual(t,[b.trainingTruth.radius],7*ones(1,c.nTrain));
verifyEqual(t,[p.forwardTruth.tau],10*ones(1,12));
verifyEqual(t,[a.truth.fStick],[p.forwardTruth.fPerm]);
verifyEqual(t,[a.truth.Dstick],[p.forwardTruth.Dperm]);
verifyEqual(t,b.regressors{1}.model.parts,["stick","ecs","sphere"]);
verifyEqual(t,a.protocol.delta,6*ones(size(protocol,1),1));
verifyGreaterThan(t,range([a.truth.radius]),0);
verifyEqual(t,[a.truth(a.trials.snr==50).radius],[a.truth(a.trials.snr==Inf).radius]);
verifyEqual(t,a.trials.b0Scale(a.trials.snr==Inf),ones(6,1));
verifyGreaterThan(t,std(a.trials.b0Scale(a.trials.snr==50)),0);

% Independently run the original exchanging-model simulator on the same seed.
ref = noise_config(["permeableStick","ecs","sphere"],protocol);
ref.delta=c.delta; ref.testPrior=c.testPrior; ref.snr=c.snr;
ref.nTissue=c.nTissue; ref.nRep=c.nRep; ref.nTrain=c.nTrain;
ref.mode="matched"; ref.fixed=struct('tau',10);
r = noise_polynomial(ref);
verifyEqual(t,r.truth,p.forwardTruth);
verifyEqual(t,r.observed,a.observed);

% Selecting a single SNR preserves the cohort and corresponding noise draws.
one = noise_exchange_mismatch(c,10,50);
keep = a.trials.snr==50;
verifyEqual(t,one.simulationSNRs,c.snr);
verifyEqual(t,one.forwardTruth,p.forwardTruth(keep));
verifyEqual(t,one.nlls.trials,a.trials(keep,:));
verifyEqual(t,one.nlls.observed,a.observed(:,:,keep));
verifyEqual(t,one.polynomial.observed,one.nlls.observed);
verifyEqual(t,one.polynomial.config.snr,50);
verifyNumElements(t,one.polynomial.regressors,1);
end

function testMatchedExchangeAtTruth(t)
% At infinite SNR a matched finite-grid prediction should reproduce every shell.
[b,Delta] = ndgrid([0 1 3 8 15],[10 40 80]);
protocol = [b(:),Delta(:),24*ones(numel(b),1)];
protocol(protocol(:,1)==0,3) = 1;
c = noise_config(["permeableStick","ecs","sphere"],protocol);
c.delta=6; c.mode="free"; c.snr=[Inf 1e6]; c.nTissue=3; c.nRep=1; c.nTrain=30; c.degree=1;
c.testPrior = struct('fractionAlpha',[4 4 2],'radius',[5 9],'p2',[.55 .8]);
r = noise_exchange_population(c,10);
a = r.nlls; poly = r.polynomial;
clean = a.trials.snr==Inf; noisy = ~clean;
verifyEqual(t,a.truth,r.forwardTruth);
verifyEqual(t,a.observed,poly.observed);
verifyEqual(t,table2array(struct2table(a.start)),table2array(struct2table(a.truth)),'AbsTol',1e-12);
verifyEqual(t,[a.truth.tau],repmat(10,1,6));
verifyGreaterThan(t,range([poly.trainingTruth.tau]),0);
verifyNumElements(t,poly.regressors{1}.targets,10);
verifyEqual(t,a.trials.fibreAxis(clean,:),a.trials.fibreAxis(noisy,:));
verifyLessThan(t,max(a.trials.startSSE(clean)),1e-24);
verifyEqual(t,table2array(struct2table(a.estimate(clean))), ...
    table2array(struct2table(a.truth(clean))),'AbsTol',1e-10);
% Verify the finite-SNR start residual is entirely accounted for by noise.
dw = protocol(:,1)>0;
deltaSignal = a.observed(dw,:,noisy)-a.observed(dw,:,clean);
verifyEqual(t,a.trials.startSSE(noisy),reshape(sum(deltaSignal.^2,[1 2]),[],1),'AbsTol',1e-19);
verifyGreaterThan(t,min(a.trials.startSSE(noisy)),0);
verifyLessThan(t,a.trials.sse(noisy),a.trials.startSSE(noisy));
end

function testExchangeBiasAndPrecision(t)
[b,Delta] = ndgrid([0 1 4 10],[10 40 80]);
protocol = [b(:),Delta(:),24*ones(numel(b),1)];
protocol(protocol(:,1)==0,3) = 1;
c = noise_config(["permeableStick","ecs","sphere"],protocol);
c.mode="free"; c.snr=100; c.nTissue=3; c.nRep=3; c.nTrain=30; c.degree=1;
c.testPrior.fractionAlpha=[4 4 2]; fixed=struct('DePar',1.6,'DePerp',.525);
testTau=[5 40]; c.trainPrior.tau=testTau;
[fits,M] = compare_exchange_assumptions(c,testTau,fixed);
ref = fits.permeable_free;
verifyEqual(t,height(M),76);  % 2 methods × (10 + 9 + 10 + 9)
for key=string(fieldnames(fits))'
    r = fits.(key);
    verifyEqual(t,r.forwardTruth,ref.forwardTruth);
    verifyEqual(t,r.exchangeTime,testTau);
    tau = [r.forwardTruth.tau];
    verifyTrue(t,all(tau>=5 & tau<=40));
    verifyGreaterThan(t,range(tau),0);
    if isfield(r.nlls.truth,'tau')
        verifyEqual(t,[r.nlls.start.tau],tau,'AbsTol',1e-12);
        verifyTrue(t,all([r.polynomial.trainingTruth.tau]>=5 & [r.polynomial.trainingTruth.tau]<=40));
    end
    verifyEqual(t,r.nlls.observed,ref.nlls.observed);
    verifyEqual(t,r.polynomial.observed,r.nlls.observed);
    verifyEqual(t,r.nlls.trials.fibreAxis,ref.nlls.trials.fibreAxis);
    if endsWith(key,"fixed_ECS")
        for k=["DePar","DePerp"]
            verifyEqual(t,[r.nlls.start.(k)],repmat(fixed.(k),1,9));
            verifyEqual(t,[r.nlls.estimate.(k)],repmat(fixed.(k),1,9));
            verifyEqual(t,[r.polynomial.estimate.(k)],repmat(fixed.(k),1,9));
            verifyEqual(t,[r.polynomial.trainingTruth.(k)],repmat(fixed.(k),1,c.nTrain));
            verifyGreaterThan(t,range([r.forwardTruth.(k)]),0);
        end
    end
    row = M(M.assumption==key & M.method=="nlls" & M.parameter=="radius",:);
    err = reshape([r.nlls.estimate.radius]-[r.nlls.truth.radius],c.nRep,[]);
    sd = sqrt(mean(var(err,0,1)));
    verifyEqual(t,row.withinTissueSD,sd,'AbsTol',1e-12);
    verifyEqual(t,row.bias,mean(err,'all'),'AbsTol',1e-12);
    verifyEqual(t,row.RMSE,sqrt(mean(err.^2,'all')),'AbsTol',1e-12);
    baseline = M(M.assumption=="permeable_free" & M.method=="nlls" & M.parameter=="radius",:);
    verifyEqual(t,row.biasShiftOverMatchedSD,(row.bias-baseline.bias)/baseline.withinTissueSD,'AbsTol',1e-10);
end
verifyTrue(t,all(isnan(M.biasOverRepeatSD(M.fixed))));
verifyEqual(t,M.withinTissueSD(M.fixed),zeros(nnz(M.fixed),1),'AbsTol',1e-14);
verifyFalse(t,any(M.parameter=="tau" & startsWith(M.assumption,"impermeable")));
end

function testIsotropicExchangeConstraints(t)
[b,Delta]=ndgrid([0 1 5 15],[10 40 80]);
protocol=[b(:),Delta(:),30*ones(numel(b),1)];
protocol(protocol(:,1)==0,3)=1;
c=noise_config(["permeableStick","ecs","sphere"],protocol,"isotropic");
c.delta=6; c.mode="free"; c.snr=50; c.nTissue=3; c.nRep=3; c.nTrain=30; c.degree=1;
c.prior.tau=[5 30]; c.prior.Dsoma=2; c.prior.p2=[.1 .5];
c.testPrior=struct('De',[.5 1.5],'Dsoma',2,'p2',[.1 .5]);
c.trainPrior=c.testPrior; c.trainPrior.tau=[5 30];
[fits,M]=compare_sandix_exchange(c,[5 30],[10 20 1e6]);
ref=fits.sandix_free; verifyEqual(t,height(M),72);
for name=string(fieldnames(fits))'
    r=fits.(name); p=r.nlls.config.prior; model=TMImix([1 20 6],c.parts,ECS="isotropic");
    fixed=isscalar(p.tau);
    verifyEqual(t,model.fitPlan(p).dimension,7-double(fixed));
    verifyEqual(t,r.forwardTruth,ref.forwardTruth);
    verifyEqual(t,r.nlls.observed,ref.nlls.observed);
    verifyEqual(t,r.polynomial.observed,r.nlls.observed);
    verifyGreaterThan(t,range([r.forwardTruth.tau]),0);
    verifyGreaterThan(t,range([r.forwardTruth.De]),0);
    verifyTrue(t,all([r.forwardTruth.tau]>=5 & [r.forwardTruth.tau]<=30));
    verifyEqual(t,[r.forwardTruth.Dsoma],2*ones(1,9));
    verifyTrue(t,all([r.forwardTruth.p2]>=.1 & [r.forwardTruth.p2]<=.5));
    for method=["nlls","polynomial"]
        verifyEqual(t,[r.(method).estimate.Dsoma],2*ones(1,9));
        verifyTrue(t,all([r.(method).estimate.p2]>=.1 & [r.(method).estimate.p2]<=.5));
    end
    verifyEqual(t,[r.polynomial.trainingTruth.Dsoma],2*ones(1,c.nTrain));
    verifyTrue(t,all([r.polynomial.trainingTruth.p2]>=.1 & [r.polynomial.trainingTruth.p2]<=.5));
    dw=r.nlls.protocol.b>0; verifyEqual(t,unique(r.nlls.protocol.Ndirs(dw)),30);
    verifyEqual(t,unique(r.nlls.protocol.L(dw)),6);
    startTruth=r.nlls.truth;
    if fixed
        [startTruth.tau]=deal(p.tau);
        verifyEqual(t,[r.polynomial.trainingTruth.tau],repmat(p.tau,1,c.nTrain));
        verifyEqual(t,[r.nlls.estimate.tau],repmat(p.tau,1,9));
        verifyEqual(t,[r.polynomial.estimate.tau],repmat(p.tau,1,9));
    end
    verifyEqual(t,table2array(struct2table(r.nlls.start)), ...
        table2array(struct2table(startTruth)),'AbsTol',1e-12);
end
verifyEqual(t,M.withinTissueSD(M.fixed),zeros(nnz(M.fixed),1));
verifyEqual(t,M.rmseRatio(M.assumption=="sandix_free" & ~M.fixed),ones(16,1));
verifyTrue(t,all(isnan(M.rmseRatio(M.parameter=="Dsoma"))));
verifyTrue(t,all(isnan(M.biasOverRepeatSD(M.fixed))));
% Exact recovery exercises the finite-direction isotropic acquisition path.
c.snr=Inf; c.nTissue=2; c.nRep=1; c.nTrain=15;
r=noise_exchange_population(c,10);
verifyLessThan(t,max(r.nlls.trials.startSSE),1e-24);
verifyEqual(t,table2array(struct2table(r.nlls.estimate)), ...
    table2array(struct2table(r.nlls.truth)),'AbsTol',1e-9);
end
