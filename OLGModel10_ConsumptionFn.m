function c=OLGModel10_ConsumptionFn(h,aprime,a,z,e,agej,Jr,r,gamma_i,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq)
% Note: these lines are essentially just a copy of the relevant part of the return fn

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
if agej<Jr
    Income=w*kappa_j*exp(gamma_i+z+e)*h+r*a; % Income is labor income and capital income
else
    Income=r*a;
end
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.

if agej<Jr % If working age
    c=(1+r)*a+(1-tau)*w*kappa_j*exp(gamma_i+z+e)*h-IncomeTax+(1+r)*AccidentBeq-aprime; % Use (z1+z2)
else % Retirement
    c=(1+r)*a+pension+(1+r)*AccidentBeq-aprime;
end


end
