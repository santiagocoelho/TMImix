function S=sphere(protocol,r,D)
% Reflecting sphere, Gaussian phase approximation, rectangular PGSE.
% z are positive roots of the derivative of spherical Bessel j_1.
persistent z
if isempty(z)
    g=@(u) (u.^2-2).*sin(u)+2*u.*cos(u);
    z=arrayfun(@(k) fzero(g,[(k-0.5)*pi,k*pi]),1:64);
end
b=protocol(:,1); big=protocol(:,2); small=protocol(:,3);
lambda=(D/r^2)*z.^2; x=small*lambda; y=(big-small)*lambda;
% Dimensionless timing factor divided by (lambda*delta)^2.
% The series also gives the analytic narrow-pulse limit at delta=0.
a=ones(size(x)); h=a; j=x>1e-3;
a(j)=-expm1(-x(j))./x(j);
h(j)=2*(x(j)+expm1(-x(j)))./x(j).^2-a(j).^2;
t=x(~j); a(~j)=1-t/2+t.^2/6-t.^3/24+t.^4/120;
h(~j)=t.*(2/3+t.*(-1/2+t.*(7/30-t/12)));
timing=h-expm1(-y).*a.^2;
adc=(2*r^2./(big-small/3)).*sum(timing./(z.^2.*(z.^2-2)),2);
S=exp(-b.*adc);
end
