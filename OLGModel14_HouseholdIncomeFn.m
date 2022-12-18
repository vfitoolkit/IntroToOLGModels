function income=OLGModel14_HouseholdIncomeFn(h,sprime,s,z,e,agej,Jr,pension,r,w,P0,D,kappa_j,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale)
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

if agej<Jr % If working age
    %consumption = labor income + accidental bequest + share holdings (including dividend) - capital gains tax - next period share holdings
    % income just is consumption but without subtracting the term for next period share holdings (-P*sprime)
    income=(1-tau_l)*w*kappa_j*exp(z+e)*Lhscale*h+((1-tau_d)*D+P0)*(s+AccidentBeq) -tau_cg*(P0-Plag)*(s+AccidentBeq); 
else % Retirement
    income=pension+((1-tau_d)*D+P0)*(s+AccidentBeq) -tau_cg*(P0-Plag)*(s+AccidentBeq);
end

% Notice that sprime>=0 is being implicitly imposed by grid on s

end
