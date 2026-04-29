function income=OLGModel14_HouseholdIncomeFn(h,sprime,s,z,e,agej,Jr,pension,w,P0,D,kappa_j,AccidentBeq,r,tau_l,tau_d,tau_cg,S_agej_first,S_agej_peak_first,S_agej_peak_last,S_agej_last)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

if sprime>=s
    agej_bought=agej;
    Plag=P0; % We are holding or buying, so no capital gains
    cg=0;
else
    if agej<=S_agej_peak_first
        agej_bought=agej-1;
        Plag=P0*(1-2*r); % Dispose of shares presumably acquired recently
    else
        % Estimate where we are past peak accumulation and mirror around to
        % proportional acquisition point
        agej_selling_pct=(agej-S_agej_peak_last)/(S_agej_last-S_agej_peak_last);
        agej_bought=S_agej_peak_first-ceil(agej_selling_pct*(S_agej_peak_first-S_agej_first));
        Plag=P0*(1-2*r)^(agej-agej_bought);
    end
    cg=tau_cg*(P0-Plag)*(s+AccidentBeq-sprime);
end

if agej<Jr % If working age
    %consumption = labor income + accidental bequest + share holdings (including dividend) - capital gains tax - next period share holdings
    % income just is consumption but without subtracting the term for next period share holdings (-P*sprime)
    income=(1-tau_l)*w*kappa_j*exp(z+e)*h+((1-tau_d)*D+P0)*(s+AccidentBeq) -cg; 
else % Retirement
    income=pension+((1-tau_d)*D+P0)*(s+AccidentBeq) -cg;
end

% Notice that sprime>=0 is being implicitly imposed by grid on s

end
