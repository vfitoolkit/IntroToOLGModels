function disposableincome=OLGModel7_DisposableIncomeFn(h,aprime,a,z,e,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A,eta1,eta2)
% Note: these lines are essentially just a copy of the relevant part of the return fn

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
if agej<Jr
    Income=w*kappa_j*exp(z+e)*h+r*a; % Income is labor income and capital income
else
    Income=r*a;
end
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.


if agej<Jr
   disposableincome=r*a+(1-tau)*w*kappa_j*exp(z+e)*h-IncomeTax; % Working age
else
   disposableincome=r*a+pension; % Retired
end
% Note: You might want +(1+r)*AccidentBeq in the disposable income
% definition, I have left it out

end
