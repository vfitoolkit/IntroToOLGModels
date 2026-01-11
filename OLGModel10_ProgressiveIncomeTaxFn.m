function IncomeTax=OLGModel10_ProgressiveIncomeTaxFn(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A,gamma_i,agej,Jr)

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
if agej<Jr % Income is labor income and capital income
    Income=w*kappa_j*exp(gamma_i+z+e)*h+r*a;
else
    Income=r*a;
end
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.

end