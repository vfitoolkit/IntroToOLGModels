function c=OLGModel4_ConsumptionFn(h,aprime,a,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A)
% Note: these lines are essentially just a copy of the relevant part of the return fn

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

if agej<Jr
    c=(1+r)*a+(1-tau)*w*kappa_j*h-aprime; % Working age
else
    c=(1+r)*a+pension-aprime; % Retired
end

end
