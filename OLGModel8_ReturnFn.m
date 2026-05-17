function F=OLGModel8_ReturnFn(h,aprime,a,z,e,sigma1,sigma2,agej,Jr,J,pension,r,A,delta,alpha,kappa_j,wg1,wg2,wg3,tau,g)
% Add deterministic growth g (compared to OLGModel6_ReturnFn)

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

F=-Inf;
if agej<Jr % If working age
    c=(1-tau)*w*kappa_j*exp(z+e)*h+(1+r)*a-(1+g)*aprime;
else % Retirement
    c=pension+(1+r)*a-(1+g)*aprime;
end

if c>0
    F=( ((c^sigma1)*((1-h)^sigma2))^(1-sigma2) )/(1-sigma2); % The utility function
end

% Warm-glow bequest
if agej==J % Final period
    warmglow=wg1*((1+(1+g)*aprime/wg2)^(1-wg3))/(1-wg3);
    F=F+warmglow;
end

end
