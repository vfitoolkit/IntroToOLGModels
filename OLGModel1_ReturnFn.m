function F=OLGModel1_ReturnFn(h,aprime,a,agej,w,sigma,psi,eta,Jr,pension,tau)

F=-Inf;

if agej<Jr
    c=(1-tau)*w*h; % Working age
else
    c=pension; % Retired
end

if c>0
    F=(c^(1-sigma))/(1-sigma) -psi*(h^(1+eta))/(1+eta);
end

end
