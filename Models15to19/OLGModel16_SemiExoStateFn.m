function prob=OLGModel16_SemiExoStateFn(emp,empprime,searcheffort,m,theta,gamma,p_jobsep)

prob=Inf; % placeholder, that will cause problems if not overwritten

% Note: job-finding probability is m*theta^(1-gamma)
p_jobfind=searcheffort*m*theta^(1-gamma);

if emp==0
    % currently non-employed
    if empprime==0
        prob=1-p_jobfind;
    elseif empprime==1
        prob=p_jobfind;
    end
elseif emp==1
    % currently employed
    if empprime==0
        prob=p_jobsep;
    elseif empprime==1
        prob=1-p_jobsep;
    end
end

end