function F=OLGModel15_ReturnFn(aprime,a,z,emp,r,h,wagemarkdown,sigma,tau,pension,lumpsum,Beq,agej,Jr)

F=-Inf;

w=h-wagemarkdown;

if agej<Jr
    c=(1+r)*a+(1-tau)*w*z*emp+lumpsum+Beq-aprime;
else
    c=(1+r)*a+pension+lumpsum-aprime;
end

if c>0
    F=(c^(1-sigma))/(1-sigma);
end




end