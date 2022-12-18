function F=OLGModel14_HouseholdReturnFn(h,sprime,s,z,e,sigma,psi,eta,agej,Jr,J,pension,r,w,P0,D,kappa_j,warmglow1,warmglow2,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

Plag=P; % As stationary general eqm

F=-Inf;
if agej<Jr % If working age
    %consumption = labor income + accidental bequest + share holdings (including dividend) - capital gains tax - next period share holdings
    c=(1-tau_l)*w*kappa_j*exp(z+e)*Lhscale*h+((1-tau_d)*D+P0)*(s+AccidentBeq) -tau_cg*(P0-Plag)*(s+AccidentBeq)-P*sprime; 
else % Retirement
    c=pension+((1-tau_d)*D+P0)*(s+AccidentBeq) -tau_cg*(P0-Plag)*(s+AccidentBeq) - P*sprime;
end

if c>0
    F=(c^(1-sigma))/(1-sigma) -psi*(h^(1+eta))/(1+eta); % The utility function
end

% Warm-glow bequest
if agej==J % Final period
    warmglow=warmglow1*(sprime^(1-warmglow2))/(1-warmglow2);
    F=F+warmglow;
end

% Notice that sprime>=0 is being implicitly imposed by grid on s

end
