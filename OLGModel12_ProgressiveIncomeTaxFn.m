function IncomeTax=OLGModel12_ProgressiveIncomeTaxFn(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2)

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
if agej<Jr
    Income=w*kappa_j1*exp(z1*e1)*h1+w*kappa_j2*exp(z2+e2)*h2+r*a; % Income is labor income and capital income
else
    Income=r*a;
end
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.

end