function F=OLGModel14_HouseholdReturnFn(h,sprime,s,z,e,sigma,psi,eta,agej,Jr,J,pension,w,P0,D,kappa_j,warmglow1,warmglow2,AccidentBeq,r,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak,S_agej_last)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize


% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

% But in fact the price does not meaningfully represent the acquisition
% cost by a younger generation that is now older in this stationary
% distribution.  However, we can use the history of acquisition and
% disposals to impute when agents are buying and selling, and thus what
% capital gains they should pay.  We imagine that stocks earn 2x the
% risk-free rate of return (i.e., 2*r) and that if we are selling before
% they peak, we are selling recently acquired stocks, whereas if we are
% selling at or after the peak of acquisition, we are selling long-term
% gains in a LIFO fashion.

% We take P0 as the price of the current stationary distribution, and we
% back-calculate what the price Plag may have been in the past.

P=P0;
if sprime>=s
    Plag=P0; % We are holding or buying, so no capital gains
else
    if agej<=S_agej_peak
        Plag=P0*(1-2*r); % Dispose of shares presumably acquired recently
    else
        % Estimate where we are past peak accumulation and mirror around to
        % proportional acquisition point
        agej_selling_pct=(agej-S_agej_peak)/(S_agej_last-S_agej_peak);
        agej_bought=S_agej_peak-ceil(agej_selling_pct*(S_agej_peak-S_agej_first));
        Plag=P0*(1-2*r)^(agej-agej_bought);
    end
end

F=-Inf;
if agej<Jr % If working age
    %consumption = labor income + accidental bequest + share holdings (including dividend) - capital gains tax - next period share holdings
    c=(1-tau_l)*w*kappa_j*exp(z+e)*h+((1-tau_d)*D+P0)*(s+AccidentBeq) -tau_cg*(P0-Plag)*(s+AccidentBeq)-P*sprime; 
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
