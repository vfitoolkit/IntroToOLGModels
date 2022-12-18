function F=OLGModel6_ReturnFn(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,J,pension,r,A,delta,alpha,kappa_j,warmglow1,warmglow2,AccidentBeq, eta1,eta2,tau)
% Add idiosyncratic shocks z and e (compared to OLGModel5_ReturnFn)

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

% Warm-glow bequest
if agej==J % Final period
    warmglow=warmglow1*(aprime^(1-warmglow2))/(1-warmglow2);
    F=F+warmglow;
end
% Notice that we have modelled the warm-glow in such a way that you only
% get it if you die in the last period (j=J). But we know there is a 1-sj
% risk of dying every period. So we might prefer to model that we get the
% warm glow bequest if we die at any age. The following commented out two lines
% implement this alternative. [note: need to add sj to inputs of ReturnFn to use it]
% warmglow=warmglowparam1*(aprime^(1-warmglowparam2))/(1-warmglowparam2); % Note: same formula as above
% F=F+(1-sj)*warmglow
% Note: if using this, have to make sure sj=0 for j=J.
% Comment: I am not aware of any study saying which of these two
% alternatives works better empirically.

end
