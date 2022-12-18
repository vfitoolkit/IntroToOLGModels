function F=OLGModel4_ReturnFn(h,aprime,a,agej,r,A,delta,alpha,sigma,psi,eta,Jr,pension,tau,kappa_j,J,warmglowparam1,warmglowparam2,AccidentBeq)
% Compared to OLGModel3_ReturnFn we have added assets, so aprime and a now
% are relevant. Comments below explain the other main changes as part of
% including assets.

% Note: in the final period (j=J) we add a warm glow bequest. 
% Without a bequest motive agents want to die with zero assets (not able to
% do so perfectly because surivial is uncertain due to sj, but will be close)
% A warm-glow bequest allows model to fit the empirical fact that people do
% still have assets when they die (and not just people who have children to
% leave assets to, people without children also have assets when they die)

% Because model has competitive labor and capital markets, and because the
% production function is Cobb-Douglas, there is a very tight relationship
% between r and w. As a result we solve the model in terms of just one of
% the two, traditionally r, and can simply get w as follows.
% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
% From K/Y, substituting the production function gives
KdivY=(KdivL^(1-alpha))/A; % Don't use this, just including it for completeness. 
% Notice that discount factor beta plays an indirect role in determining K/Y via determining r
% We know w=MPL (MPL is marginal product of labour)
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

F=-Inf;

if agej<Jr
    c=(1+r)*a+(1-tau)*w*kappa_j*h+(1+r)*AccidentBeq-aprime; % Working age
else
    c=(1+r)*a+pension+(1+r)*AccidentBeq-aprime; % Retired
end

if c>0
    F=(c^(1-sigma))/(1-sigma) -psi*(h^(1+eta))/(1+eta);
end

% Warm-glow bequest
if agej==J % Final period
    warmglow=warmglowparam1*(aprime^(1-warmglowparam2))/(1-warmglowparam2);
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
