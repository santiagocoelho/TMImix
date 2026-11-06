function tests=test_physics
tests=functiontests(localfunctions);
end

function setupOnce(t)
t.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));
end

function testSMIProjection(t)
b=[0;0.7;2;6]; m=TMImix([b,24+0*b,4+0*b],["stick","ecs","water"]);
p=m.defaults(); p.fStick=.4; p.fECS=.5; p.fWater=.1; p.p2=.63;
S=m.invariants(p);
for j=1:numel(b)
    k=@(x) .4*exp(-b(j)*p.Dstick*x.^2)+.5*exp(-b(j)*(p.DePerp+ ...
        (p.DePar-p.DePerp)*x.^2))+.1*exp(-b(j)*p.Dwater);
    verifyEqual(t,S(j,:),[integral(k,0,1),p.p2*abs(integral(@(x) k(x).*(3*x.^2-1)/2,0,1))], ...
        'AbsTol',1e-12);
end
end

function testSphere(t)
q=[1 18 4;2 30 8;4 55 12]; m=TMImix(q,"sphere"); p=m.defaults();
% Independent fixtures from SANDIX fun_comps_ave, 31 eigenmodes, R=7, D0=2.
verifyEqual(t,m.invariants(p),[[.642648542646893;.629317701697557;.657000626112385],zeros(3,1)],'AbsTol',5e-10);
A=tmix.sphere([2 30 0;2 30 1e-7],7,2);
verifyEqual(t,A(1),A(2),'AbsTol',1e-8);
% Fully restricted and free diffusion limits.
verifyEqual(t,tmix.sphere([2 20 5],.05,2),1,'AbsTol',1e-6);
verifyGreaterThan(t,tmix.sphere([2 20 5],100,2),exp(-4));
end

function testExchangeAgainstODE(t)
q=[3 28 9]; m=TMImix(q,["ecs","permeableStick"],Steps=32);
p=m.defaults(); p.fPerm=.37; p.fECS=.63; p.tau=13;
u=[.1 .55 .9]; ref=zeros(size(u)); eta=p.fPerm;
R=[-(1-eta) eta;1-eta -eta]/p.tau;
for j=1:numel(u)
    D=diag([p.Dperm*u(j)^2,p.DePerp+(p.DePar-p.DePerp)*u(j)^2]);
    q2=@(s) q(1)/(q(2)-q(3)/3)*(max(0,min([1,s/q(3),(q(2)+q(3)-s)/q(3)])))^2;
    [~,Y]=ode45(@(s,y) (R-q2(s)*D)*y,[0 q(3) q(2) q(2)+q(3)], ...
        [eta;1-eta],odeset('RelTol',1e-11,'AbsTol',1e-13));
    ref(j)=sum(Y(end,:));
end
verifyEqual(t,m.kernel(p,u),ref,'AbsTol',2e-5);
mn=TMImix(q,m.parts,Exchange="narrow"); E=mn.kernel(p,u);
for j=1:numel(u)
    D=diag([p.Dperm*u(j)^2,p.DePerp+(p.DePar-p.DePerp)*u(j)^2]);
    verifyEqual(t,E(j),sum(expm((q(2)-q(3)/3)*R-q(1)*D)*[eta;1-eta]),'AbsTol',1e-13);
end
verifyGreaterThan(t,max(abs(ref-E)),1e-4);
end

function testExchangeLimitsAndAllParts(t)
q=[0 25 4;2 20 4;2 55 4]; u=[0 .4 1]; m=TMImix(q,["water","dot","sphere","stick","ecs","permeableStick"]);
p=m.defaults(); p.tau=Inf;
% Reconstruct a no-exchange mixture from individual kernels.
expected=zeros(3,3);
for part=["water","dot","sphere","stick","ecs"]
    solo=TMImix(q,part); a=solo.defaults();
    for k=string(fieldnames(a))', a.(k)=p.(k); end
    a.(solo.fractions)=1; expected=expected+p.(solo.fractions)*solo.kernel(a,u);
end
expected=expected+p.fPerm*exp(-q(:,1)*(p.Dperm*u.^2));
verifyEqual(t,m.kernel(p,u),expected,'AbsTol',1e-12);
p.tau=20; verifyEqual(t,m.kernel(p,u)*[1;0;0],m.kernel(p,0),'AbsTol',1e-12);
m=TMImix(q,["permeableStick","ecs"]); p=m.defaults(); S=m.invariants(p);
verifyLessThan(t,S(3,1),S(2,1)); p.fECS=0; p.fPerm=1;
verifyEqual(t,m.kernel(p,u),exp(-q(:,1)*(p.Dperm*u.^2)),'AbsTol',1e-13);
end

function testFODFAndMeasuredInvariants(t)
q=[.5 16 4;2 40 4]; m=TMImix(q,["stick","ecs","sphere"]);
p=m.defaults(); d=fibonacci(250); a=[1 0 0; .5 sqrt(.75) 0; 0 0 1]; w=[.5;.4;.1];
p.p2=TMImix.odfP2(a,w);
verifyEqual(t,TMImix.odfP2(eye(3),[1;1;1]),0,'AbsTol',1e-14);
Y=m.signal(p,d,a,w);
R=[0 -1 0; 1 0 0; 0 0 1];
verifyEqual(t,m.signal(p,d*R,-a*R,2*w),Y,'AbsTol',1e-14);
W=tmix.projector(d,12); verifyEqual(t,tmix.measure(Y,W),m.invariants(p),'AbsTol',2e-6);
% Convolution also works across node/gradient block boundaries.
verifyEqual(t,m.signal(p,d(1:35,:),repmat(a,45,1),repmat(w,45,1)),Y(1:35,:),'AbsTol',1e-13);
end

function d=fibonacci(n)
i=(0:n-1)'; z=(i+.5)/n; phi=pi*(3-sqrt(5))*i;
d=[sqrt(1-z.^2).*cos(phi),sqrt(1-z.^2).*sin(phi),z];
end
