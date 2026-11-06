function W=projector(dirs,L)
% Rows map directional DWI onto the first six real SH coefficients.
% Fit nuisance orders through L to reduce leakage into l=0,2; discard them.
if nargin<2, L=4; end
tmix.unitRows(dirs);
validateattributes(L,{'numeric'},{'scalar','integer','>=',2});
assert(mod(L,2)==0,'L must be even.');
phi=atan2(dirs(:,2),dirs(:,1)); u=max(-1,min(1,dirs(:,3)))';
Y=zeros(numel(u),(L+1)*(L+2)/2); j=0;
for l=0:2:L
    P=legendre(l,u);
    for m=-l:l
        a=abs(m); j=j+1;
        v=sqrt((2*l+1)/(4*pi)*exp(gammaln(l-a+1)-gammaln(l+a+1)))*P(a+1,:)';
        if m<0, v=v.*(sqrt(2)*sin(a*phi));
        elseif m>0, v=v.*(sqrt(2)*cos(a*phi)); end
        Y(:,j)=v;
    end
end
assert(rank(Y)==size(Y,2),'Too few independent directions for this SH order.');
W=pinv(Y); W=W(1:6,:);
end
