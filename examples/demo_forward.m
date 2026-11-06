% All six compartments, two crossing fibre bundles, two pulse timings.
m=TMImix([0 20 5; 2 20 5; 2 50 5], ...
    ["stick","ecs","permeableStick","sphere","dot","water"]);
p=m.defaults();
p.fStick=.23; p.fECS=.35; p.fPerm=.20; p.fSphere=.15; p.fDot=.02; p.fWater=.05;
axes=[1 0 0; cosd(65) sind(65) 0]; mass=[0.65;0.35];
p.p2=TMImix.odfP2(axes,mass);
directional=m.signal(p,eye(3),axes,mass)  % directions x shells
invariants=m.invariants(p)                % shells x [S0,S2]
