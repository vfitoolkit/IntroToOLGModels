function IncomeTax=OLGModel5_ProgressiveIncomeTaxFn(h,aprime,a,eta1,eta2,kappa_j,r,delta,alpha,A)

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
Income=w*kappa_j*h+r*a; % Income is labor income and capital income
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.

end