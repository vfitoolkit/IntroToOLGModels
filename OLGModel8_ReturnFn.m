function F=OLGModel8_ReturnFn(h,aprime,a,z,e,sigma1,sigma2,agej,Jr,pension,r,A,delta,alpha,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj,tau,g)
% Add deterministic growth g (compared to OLGModel6_ReturnFn)

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

F=-Inf;
if agej<Jr % If working age
    c=(1-tau)*w*kappa_j*exp(z+e)*h+(1+r)*a-(1+g)*aprime;
else % Retirement
    c=pension+(1+r)*a-(1+g)*aprime;
end

if c>0
    F=( ((c^sigma1)*((1-h)^sigma2))^(1-sigma2) )/(1-sigma2); % The utility function
end

% add the warm glow to the return, but only near end of life
if agej>=Jr+10
    % Warm glow of bequests
    warmglow=warmglow1*(((1+g)*aprime-warmglow2)^(1-warmglow3))/(1-warmglow3);
    % Modify for beta and sj (get the warm glow next period if die)
    warmglow=beta*(1-sj)*warmglow;
    % add the warm glow to the return
    F=F+warmglow;
end

end
