function c=OLGModel12_ConsumptionFn(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,pension,r,kappa_j1,kappa_j2,gamma_1,gamma_2,A,alpha,delta,eta1,eta2,tau,AccidentBeq)
% Note: these lines are essentially just a copy of the relevant part of the return fn

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
if agej<Jr % Income is labor income and capital income
    Income=w*kappa_j1*exp(gamma_1+z1+e1)*h1+w*kappa_j2*exp(gamma_2+z2+e2)*h2+r*a;
else
    Income=r*a;
end
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.

if agej<Jr % If working age use (z1+z2)
    c=(1+r)*a+(1-tau)*(w*kappa_j1*exp(gamma_1+z1+e1)*h1+w*kappa_j2*exp(gamma_2+z2+e2)*h2)-IncomeTax+(1+r)*AccidentBeq-aprime;
else % Retirement
    c=(1+r)*a+pension+(1+r)*AccidentBeq-aprime;
end

end
