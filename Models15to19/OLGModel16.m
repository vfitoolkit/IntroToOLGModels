%% OLG Model 16: Endogenous Search Effort
% Mostly same as OLG Model 15, but endogenising the search effort.
% Main change this makes for the code is that we need search effort as a
% decision variable, and that the employed/non-employed state is now as
% semi-exogenous state (previously it was an exogenous markov).

%% Action and State Space
n_d=3; % search intensity (hardcoded n_d=3 into the ReturnFn and d_grid)
n_a=1001; % assets
n_semiz=2; % emp, not-employed/employed (0/1)
n_z=31; % z, idiosyncratic labor productivity
N_j=81; % ages 20 to 100

% Note: to get GeneralEqmEqns near zero, you need a variety of people
% making different choices on search effort. In this model that requires
% that you have a large number of points on z. [Adding to n_d only helps 
% very little, as it is about the expected value of the firm which is loosely
% an expectation of the z over distribution of unemployed, and the d only 
% changes the probabilities a little.]

%% Parameters

% Age
Params.agej=1:1:N_j;
Params.agejshifter=19; % first model period is age 20
Params.Jr=65-Params.agejshifter; % Retirement age

% Discount factor
Params.beta=0.95;

% Preferences
Params.sigma=2; % CES utility parameter (risk aversion)
Params.psi1=0.02; % disutility of searcheffort>0
Params.psi2=0.7; % curvature of disutility of searcheffort

% Prices
Params.r=0.05; % interest rate, initial guess, determined in general eqm
Params.h=1.2; % marginal product of labor, initial guess, determined in general eqm
Params.wagemarkdown=0.5; % marginal product of labor - wage; initial guess, determined in general eqm
% wage: w=h-wagemarkdown;

% Taxes
Params.tau=0.06;

% Pensions
Params.pension=0.1;

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report, volume 58, number 10, March 2010."
% I took them from first column (qx) of Table 1 (Total Population)
% Conditional death probabilities
Params.dj=[0.006879, 0.000463, 0.000307, 0.000220, 0.000184, 0.000172, 0.000160, 0.000149, 0.000133, 0.000114, 0.000100, 0.000105, 0.000143, 0.000221, 0.000329, 0.000449, 0.000563, 0.000667, 0.000753, 0.000823,...
    0.000894, 0.000962, 0.001005, 0.001016, 0.001003, 0.000983, 0.000967, 0.000960, 0.000970, 0.000994, 0.001027, 0.001065, 0.001115, 0.001154, 0.001209, 0.001271, 0.001351, 0.001460, 0.001603, 0.001769, 0.001943, 0.002120, 0.002311, 0.002520, 0.002747, 0.002989, 0.003242, 0.003512, 0.003803, 0.004118, 0.004464, 0.004837, 0.005217, 0.005591, 0.005963, 0.006346, 0.006768, 0.007261, 0.007866, 0.008596, 0.009473, 0.010450, 0.011456, 0.012407, 0.013320, 0.014299, 0.015323,...
    0.016558, 0.018029, 0.019723, 0.021607, 0.023723, 0.026143, 0.028892, 0.031988, 0.035476, 0.039238, 0.043382, 0.047941, 0.052953, 0.058457, 0.064494,...
    0.071107, 0.078342, 0.086244, 0.094861, 0.104242, 0.114432, 0.125479, 0.137427, 0.150317, 0.164187, 0.179066, 0.194979, 0.211941, 0.229957, 0.249020, 0.269112, 0.290198, 0.312231, 1.000000]; 
% dj covers Ages 0 to 100
Params.sj=1-Params.dj(21:101); % Conditional survival probabilities
Params.sj(end)=0; % We did not need this in the previous model, but now it is important to help keep track of accidental bequests (otherwise we would have to use a more complicated formula)


% Labor market
% Households and intermediate firms interact in search-and-matching labor market
% Matching function: M(U,V)=m U^gamma V^(1-gamma)
% Market tightness: theta=V/U.
Params.m=1; % Matching fn technology level
Params.theta=0.8; % market tightness (ratio of vacancies to not-employed)
                   % initial guess, theta is determined in general eqm
Params.gamma=0.7; % exponent in the matching function

% Labor hire
% Intermediate firms rent labor to final goods firm
Params.eta=0.5; % workers wage bargaining power
Params.c=0.2; % cost of posting vacancies
Params.omega=0.1; % workers outside option

% Lump sum the profits of intermediate firms to households
Params.lumpsum=0.05;
% Accidental bequests
Params.Beq=0.03;

% Idiosyncratic labor producitivity
Params.rho_z=0.7; % autocorrelation
Params.sigma_z_epsilon=0.1; % std dev of innovations to z

% Production fn
Params.alpha=0.3;
Params.delta=0.05;

%% Grids

d_grid=linspace(0,1,n_d)'; % search effort [because of how we set the SemiExoStateFn, d must be 0 to 1]

% asset grid
a_grid=3*linspace(0,1,n_a)'.^3; %.^3 puts more points near zero

% To be able to determine the transition parameters of the exogneous markov
% in general eqm we have to set it up using ExogShockFn
% The idiosyncratic labor productivity does not need to be determined in
% general eqm, but we have to do both shocks together
% Params.p_jobfind=Params.m*Params.theta^(1-Params.gamma); % job finding probability
Params.p_jobsep=0.2; % job separation probability

% idiosyncratic worker productivity
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,Params.rho_z,Params.sigma_z_epsilon,n_z);
z_grid=exp(z_grid);

%% Semi-exogenous state
semiz_grid=[0; 1]; % 0=non-employed, 1=employed

% Define the transition probabilities of the semi-exogenous states
vfoptions.SemiExoStateFn=@(emp,empprime,searcheffort,m,theta,gamma,p_jobsep) OLGModel16_SemiExoStateFn(emp,empprime,searcheffort,m,theta,gamma,p_jobsep);
% It is hardcoded that only the 'last' decision variable can influence the transition probabilities of the semi-exogenous states
% The semi-exogenous states must be included in the return fn, fns to evaluate, etc. The ordering must be that semi-exogenous states come after this period endogenous state(s) and before any markov exogenous states, so (...,a,semiz,z,...)

% Put everything into vfoptions and simoptions
vfoptions.n_semiz=n_semiz;
vfoptions.semiz_grid=semiz_grid;
simoptions.n_semiz=vfoptions.n_semiz;
simoptions.semiz_grid=vfoptions.semiz_grid;
simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
% need to put d_grid in simoptions when using semi-exo state
simoptions.d_grid=d_grid;

% Note: theta must be 0 to 1, so when we get to general eqm we will constrain this.


%% Discount factor and Return fn
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(searcheffort,aprime,a,emp,z,r,h,wagemarkdown,sigma,psi1,psi2,tau,pension,lumpsum,Beq,agej,Jr)...
    OLGModel16_ReturnFn(searcheffort,aprime,a,emp,z,r,h,wagemarkdown,sigma,psi1,psi2,tau,pension,lumpsum,Beq,agej,Jr);

%% Divide-and-conquer and grid interpolation layer
vfoptions.divideandconquer=1; % turn on divide-and-conquer
vfoptions.gridinterplayer=1; % turn on grid interpolation layer
vfoptions.ngridinterp=20; % 20 evenly-spaced points between each pair of consecutive a_grid points
simoptions.gridinterplayer=vfoptions.gridinterplayer; % grid interpolation layer must also be set in simoptions
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Solve value fn and policy fn
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

%% Initial dist
jequaloneDist=zeros([n_a,n_semiz,n_z]);
jequaloneDist(1,:,3)=shiftdim([0.8,0.2],-1); % everyone is born with zero assets, 80% not-employed, and median productivity.

AgeWeightParamNames={'mewj'};
Params.mewj=ones(1,N_j)/N_j;

%% Stationary dist
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

%% FnsToEvaluate
FnsToEvaluate.L=@(searcheffort,aprime,a,emp,z) z*(emp==1); % effective labor units of the employed [note: retired all have e=0, so no need to remove them]
FnsToEvaluate.K=@(searcheffort,aprime,a,emp,z) a; % assets
FnsToEvaluate.pensionrevenue=@(searcheffort,aprime,a,emp,z,tau,h,wagemarkdown,agej,Jr) tau*(h-wagemarkdown)*z*emp*(agej<Jr); % assets
FnsToEvaluate.pensionspending=@(searcheffort,aprime,a,emp,z,pension,agej,Jr) pension*(agej>=Jr); % assets
FnsToEvaluate.BeqLeft=@(searcheffort,aprime,a,emp,z,sj) aprime*(1-sj); % assets
FnsToEvaluate.BeqReceived=@(searcheffort,aprime,a,emp,z,Beq,agej,Jr) Beq*(agej<Jr); % assets
FnsToEvaluate.FirmProfits=@(searcheffort,aprime,a,emp,z,wagemarkdown) wagemarkdown*z*emp; % assets

AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

%% Better starting guess
Params.h=0.90;
Params.r=0.06;
Params.wagemarkdown=0.28;
Params.theta=0.78;
Params.lumpsum=0.13;
Params.pension=0.12;
Params.Beq=0.03;

%% General eqm
GEPriceParamNames={'h','r','wagemarkdown','theta','lumpsum','pension','Beq'};
% Note: we could eliminate h and laborhire from GE, and instead take
% advantage of the fact that h can be expressed in terms of r (standard
% trick for stationary eqm with Cobb-Douglas prodn fn).

GeneralEqmEqns.capitalmarket=@(r,K,L,alpha,delta) r- (alpha*(K^(alpha-1))*(L^(1-alpha))-delta);
GeneralEqmEqns.laborhire=@(h,K,L,alpha) h-(1-alpha)*(K^(alpha))*(L^(-alpha));
GeneralEqmEqns.intermediatefirm=@(wagemarkdown,h,eta,omega,c,theta) (h-wagemarkdown)-((1-eta)*omega+eta*(h+c*theta));
GeneralEqmEqns.firmfreeentrycondn=@(c,m,theta,gamma,r,Efirmvalue) c/(m*theta^(-gamma))-(1/(1+r))*Efirmvalue;
GeneralEqmEqns.firmprofits=@(lumpsum,FirmProfits) lumpsum-FirmProfits; % note: because everyone gets lumpsum, we can just use it here rather than take AggVar (as it would just be intergral of a constant over a mass of one)
GeneralEqmEqns.pensionbalance=@(pensionrevenue,pensionspending) pensionrevenue-pensionspending;
GeneralEqmEqns.bequests=@(BeqLeft,BeqReceived) BeqLeft-BeqReceived;
% I think of it intuitively like:
% marginal product of capital determines r [final good producer + competitive markets]
% marginal product of labor determines h [final good producer + competitive markets]
% intermediate firms determine wagemarkdown [markup h over w (so markdown w from h)]
% firm free entry condition determines theta [entry increases V and thus theta, happens until balanced]
% The others are just about cleaning up.
% Of course, in reality all the prices determine all the general eqm eqns jointly.

heteroagentoptions.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions)...
    OLGModel16_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,heteroagentoptions,vfoptions,simoptions);
% Note: these inputs cannot be modified
% Note: any other input you want, just 'hide' a copy in heteroagentoptions

heteroagentoptions.constrain0to1={'theta'}; % otherwise the transition probabilities for employment semi-exo state are not 0 to 1

heteroagentoptions.verbose=1;
% heteroagentoptions.fminalgo=8; % lsqnonlin(), this requires Matlab Optimization Toolbox
heteroagentoptions.fminalgo=1; % fminsearch()
heteroagentoptions.toleranceGEcondns=10^(-4);
[p_eqm,~,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

%% Analyse solution
% Update to general eqm params
Params.r=p_eqm.r;
Params.h=p_eqm.h;
Params.wagemarkdown=p_eqm.wagemarkdown;
Params.theta=p_eqm.theta;
Params.lumpsum=p_eqm.lumpsum;
Params.pension=p_eqm.pension;
Params.Beq=p_eqm.Beq;
% Some things of interest
Params.w=Params.h-Params.wagemarkdown;
Params.p_jobfind=Params.m*Params.theta^(1-Params.gamma); % job finding probability
Params.p_jobfill=Params.m*Params.theta^(-Params.gamma); % job filling probability

[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
FnsToEvaluate2=FnsToEvaluate;
FnsToEvaluate2.earnings=@(searcheffort,aprime,a,emp,z,w) w*z*emp;
FnsToEvaluate2.employmentstatus=@(searcheffort,aprime,a,emp,z) emp;
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate2,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate2,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);



StationaryDistj=squeeze(sum(sum(sum(StationaryDist,1),2),3));
StationaryDistjj=reshape(StationaryDist,[n_a*n_semiz*n_z,N_j]);
min(StationaryDistjj,[],1)


