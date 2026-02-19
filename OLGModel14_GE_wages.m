function F=OLGModel14_GE_wages(L_h,L_f,w)
% A continuous function that returns zero when L_h is balanced with L_f
% and otherwise attempts to steer wages (w) toward balancing L_h and L_f

F=0;

if (L_h-L_f)^2 <= 10^(-6)
    return
elseif L_f>L_h
    if w>=1
        % Small increase
        F=(L_f/L_h)*(w-1);
    else
        % Larger increase
        F=(L_f/L_h)^2*(1-w);
    end
else
    if w>1
        % Large decrease (negative)
        F=(L_h/L_f)^2*(1-w);
    else
        % Smaller decrease (negative)
        F=(L_h/L_f)*(w-1);
    end
end

end