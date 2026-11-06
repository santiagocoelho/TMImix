function E=exchange(protocol,Di,De,eta,tau,mode,n)
% Two-site Karger evolution at detailed balance, with pulse-resolved q(t)^2.
b=protocol(:,1); big=protocol(:,2); small=protocol(:,3);
if isinf(tau)
    E=eta*exp(-b*Di)+(1-eta)*exp(-b*De); return
end
t=big-small/3; I=eta; O=1-eta;
if mode=="narrow" || all(small==0)
    [I,O]=step(I,O,b,t,Di,De,eta,tau);
else
    q2=b./t; h=small/n;
    % Exact integral of q^2 on each ramp subinterval; exact flat top.
    scale=((1:n).^2-(1:n)+1/3)/n^2;
    for s=scale, [I,O]=step(I,O,(q2.*h)*s,h,Di,De,eta,tau); end
    [I,O]=step(I,O,q2.*(big-small),big-small,Di,De,eta,tau);
    for s=fliplr(scale), [I,O]=step(I,O,(q2.*h)*s,h,Di,De,eta,tau); end
end
E=I+O;
end

function [I,O]=step(I,O,b,h,Di,De,eta,tau)
% exp(-A) for a 2x2 block; no exp(+large) or division by an eigenvalue gap.
i=b*Di; e=b*De; ki=h*((1-eta)/tau); ke=h*(eta/tau);
a=i+ki; d=e+ke;
gap=hypot(a-d,2*sqrt(ki.*ke));
fast=(a+d+gap)/2;
slow=(i.*e+i.*ke+e.*ki)./max(fast,realmin);
v=exp(-slow); c=-expm1(-gap)./max(gap,realmin); c(gap==0)=1; c=v.*c;
newI=(v+c.*(slow-a)).*I+(c.*ke).*O;
O=(c.*ki).*I+(v+c.*(slow-d)).*O;
I=newI;
end
