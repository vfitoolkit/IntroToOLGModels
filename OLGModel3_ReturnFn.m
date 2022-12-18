function F=OLGModel3_ReturnFn(h,aprime,a,agej,w,sigma,psi,eta,Jr,pension,tau,kappa_j)
% Note: VFI Toolkit automatically deals with dependence of kappa_j (and
% agej) on the age j, and the input to the return function is just the
% relevant value for the 'current' age.

F=-Inf;

if agej<Jr
    c=(1-tau)*w*kappa_j*h; % Working age
else
    c=pension; % Retired
end

if c>0
    F=(c^(1-sigma))/(1-sigma) -psi*(h^(1+eta))/(1+eta);
end

end
