function c=OLGModel5_ConsumptionFn(h,aprime,a,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A,eta1,eta2)
% Note: these lines are essentially just a copy of the relevant part of the return fn

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
Income=w*kappa_j*h+r*a; % Income is labor income and capital income
IncomeTax=eta1+eta2*log(Income)*Income;

if agej<Jr
    c=(1+r)*a+(1-tau)*w*kappa_j*h-IncomeTax-aprime; % Working age
else
    c=(1+r)*a+pension-aprime; % Retired
end

end
