function cg=OLGModel14_HouseholdCapitalGainsFn(h,sprime,s,z,e,agej,P0,AccidentBeq,r,tau_cg,S_agej_first,S_agej_peak,S_agej_last)
% Replace assets with 'share holdings'
% Get rid of progressive taxes
% Add Lhnormalize

% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

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

cg=tau_cg*(P0-Plag)*(s+AccidentBeq);

end
