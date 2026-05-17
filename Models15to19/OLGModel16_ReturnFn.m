function F=OLGModel16_ReturnFn(searcheffort,aprime,a,emp,z,r,h,wagemarkdown,sigma,psi1,psi2,tau,pension,lumpsum,Beq,agej,Jr)

F=-Inf;

w=h-wagemarkdown;

if agej<Jr
    c=(1+r)*a+(1-tau)*w*z*emp+lumpsum+Beq-aprime;
else
    c=(1+r)*a+pension+lumpsum-aprime;
end

if c>0
    F=(c^(1-sigma))/(1-sigma)-psi1*(searcheffort>=0)*searcheffort^psi2;
end




end