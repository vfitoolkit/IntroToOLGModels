function F=OLGModel9B_ReturnFn(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,pension,r,A,delta,alpha,kappa_j,AccidentBeq, eta1,eta2,tau)
% Just deletes the warm-glow of bequests (compared to OLGModel6_ReturnFn),
% as this has to be handled specially for EZ preferences.

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


F=-Inf;
if agej<Jr % If working age
    c=(1+r)*a+(1-tau)*w*kappa_j*exp(z+e)*h-IncomeTax+(1+r)*AccidentBeq-aprime; % Use (z1+z2)
else % Retirement
    c=(1+r)*a+pension+(1+r)*AccidentBeq-aprime;
end

if c>0
    F=(c^(1-sigma))/(1-sigma) -psi*(h^(1+eta))/(1+eta); % The utility function
end


end
