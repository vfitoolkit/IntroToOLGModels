function CustomStats=OLGModel16_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions)
% CustomModelStats: you CANNOT change the inputs at all.
% Note: any other input you want, just 'hide' a copy in heteroagentoptions
CustomStats=struct();
% CustomStats can contain any 'parameters' you want to create.
% It must be a structure, and the fields must be the 'parameter names'
% E.g., CustomStats.sigma
% Note, this is the same kind of setup as Params, the usual parameter struture

%% OLG Model 16
% NO REAL DIFFERENCE. JUST ABOUT FACT THAT pi_z IS NOW JUST PRODUCTIVITY,
% AND StationaryDist NOW has (...emp,z,...)

% The custom stat we need is the for the free-entry general eqm eqn.
% We need to calculate the 'expected value of entry', which is the integral
% of the firm value function, with respect to the distribution of
% not-employed households.

% This model is simple enough, that we can get a closed-form expression for
% the firm value fn.
% StationaryDist for z==0 is the distribution of not-employed households

z_grid=z_gridvals_J(:,1); % only one z,  and is independent of age
pi_z=pi_z_J(:,:,1); % only one z, and is independent of age

%% First, calculate J (intermediate firm value function)
% We could do this by having firms as a permanent type in the model.
% But since there are no 'decisions' for firm, we can just do it here
% 'manually'.

% We have
% w=h-wagemarkdown;
% F=(h-w)*z1_grid;
% But this simplifies to
F=Params.wagemarkdown*z_grid;

J=zeros(n_z(1),Params.Jr-1); % J(z,j)
J(:,Params.Jr-1)=F; % recall, job-separation probability =1 in period before retirement (because worker will retire)
for jjrev=1:Params.Jr-2
    jj=Params.Jr-1-jjrev;
    dj=1-Params.sj(jj); % probability of dying
    firmdiscountfactor=1/(1+Params.r); % firm discounts using market rate of interest
    J(:,jj)=F+(1-(Params.p_jobsep+dj-Params.p_jobsep*dj))*firmdiscountfactor*(pi_z*J(:,jj+1));
end

% Take expectation of J w.r.t. the household distribution
% J is z-by-(Jr-1)
% StationaryDist is (a,z,e,j)
EffortWeightedStationaryDist=shiftdim(d_grid(Policy(1,:,:,:,:)),1).*StationaryDist;
Distofunemployed=squeeze(sum(EffortWeightedStationaryDist(:,1,:,1:Params.Jr-1),1)); % now: (z,j)
% StationaryDist(:,1,:,:) is unemployed, and then sum the a dimensions (for working age)
% So Distofunemployed is the search-effort-weighted distribution of the unemployed, over z and j (but only working ages)

% So we want
CustomStats.Efirmvalue=sum(sum(J.*Distofunemployed));
% Done!


end