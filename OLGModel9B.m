%% OLG Model 9B: Epstein-Zin preferences
% Modify OLGModel6 to use Epstein-Zin preferences.
% Epstein-Zin preferences seperate the 'intertemporal elasticity of substitution' from the 'risk aversion' 
% (a single parameter determines both in standard vonNeumann-Morgenstern preferences).
% This example uses Epstein-Zin preferences in utility-units (rather than the more traditional
% consumption-units formulation of EZ prefs)

% There are essentially three parts to using Epstein-Zin preferences.
% 1. Use vfoptions to state that you are using Epstein-Zin preferences.
% 2. Set the appropriate preference parameters
% 3. Be careful about warm-glow of bequests (remove from return function)

% 1. Use vfoptions to state that you are using Epstein-Zin preferences.
vfoptions.exoticpreferences='EpsteinZin'; % Use Epstein-Zin preferences
vfoptions.EZpositiveutility=0; % utility function is negative

% 2. Set the appropriate preference parameters
vfoptions.EZriskaversion='phi';
Params.phi=2; % Risk aversion
% This parameter controls the 'additional risk aversion' relative to vNM risk preferences.

% 3. Be careful about warm-glow of bequests (remove from return function)
% To be able to use a warm-glow-of-bequests with Epstein-Zin preferences we have to distingush
% conditional survival probabilities from the regular discount factor. 
% So below the discount factor is now just beta, and we add
vfoptions.survivalprobability='sj';
% We no longer include warm-glow of bequests in the return fn.

% 4. Warm-glow of bequests
% Using warm-glow of bequests together with EZ preferences is subtle, so dealing with them has been mostly automated.
% Need to define two things (if you don't want bequests you simply do not define these)
vfoptions.WarmGlowBequestsFn=@(aprime,sigma,wg,agej,J) (agej==J)*wg*(aprime^(1-sigma))/(1-sigma); % First input arguement must be aprime, after than can be any parameters
Params.wg=1; % controls strength of bequest motive
% Comment: Loosely speaking you want the WarmGlowBequestsFn to output the 'same'
% thing as the return fn. Our utility function has (c^(1-sigma))/(1-sigma)
% and hence we set the warmglow to (aprime^(1-sigma))/(1-sigma). We can
% then control the importance of the warm-glow of bequests by multiplying
% it by a constant, here called wg. Note that to keep this in line with previous models we
% also include a term so that the warm-glow of bequests in only non-zero
% in the final period, hence the (agej==J)
% Comment: If a parameter in the WarmGlowBequestsFn depends on age, then it is the last period of life
% from which the parameter value is taken. So be careful. This may mean you want to create
% an offset version of a parameter to put into the WarmGlowBequestsFn. (E.g., if
% you die at the end of period 20, then it is the period 20 parameter
% values that will be used to evaluate WarmGlowBequestsFn; people often 
% think of the warm-glow as being received in the following period.)

% Have also changed the ReturnFn below (just deleting the warm-glow of bequests from it) 

% Done! Everything else below is unchanged :)

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=51; % Endogenous labour choice (fraction of time worked)
n_a=301; % Endogenous asset holdings
% Exogenous labor productivity units shocks (next two lines)
n_z=15; % AR(1) with age-dependent params
vfoptions.n_e=3; % iid
N_j=Params.J; % Number of periods in finite horizon

%% Parameters

% Discount rate
Params.beta = 0.96;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 10; % Weight on leisure

Params.A=1; % Aggregate TFP. Not actually used anywhere.
% Production function
Params.alpha = 0.3; % Share of capital
Params.delta = 0.1; % Depreciation rate of capital

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr=46;
% Population growth rate
Params.n=0.02; % percentage rate (expressed as fraction) at which population growths

% Age-dependent labor productivity units
Params.kappa_j=[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];
% Life-cycle AR(1) process z, on (log) labor productivity units
% Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Originals just cover ages 24 to 60, so I create these, and then repeat first and last periods to fill it out
Params.rho_z=0.7596+0.2039*((1:1:37)/10)-0.0535*((1:1:37)/10).^2+0.0028*((1:1:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
Params.sigma_epsilon_z=0.0518-0.0405*((1:1:37)/10)+0.0105*((1:1:37)/10).^2-0.0002*((1:1:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Note that 37 covers 24 to 60 inclusive
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
Params.rho_z=[Params.rho_z(1)*ones(1,4),Params.rho_z,Params.rho_z(end)*ones(1,4),zeros(1,100-65+1)];
Params.sigma_epsilon_z=[Params.sigma_epsilon_z(1)*ones(1,4),Params.sigma_epsilon_z,Params.sigma_epsilon_z(end)*ones(1,4),Params.sigma_epsilon_z(end)*ones(1,100-65+1)];
% Transitory iid shock
Params.sigma_e=0.0410+0.0221*((24:1:60)/10)-0.0069*((24:1:60)/10).^2+0.0008*((24:1:60)/10).^3;
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
Params.sigma_e=[Params.sigma_e(1)*ones(1,4),Params.sigma_e,Params.sigma_e(end)*ones(1,4),Params.sigma_e(end)*ones(1,100-65+1)];
% Note: These will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)
% Note: Karahan & Ozkan (2013) also have a fixed effect (which they call alpha) and which I ignore here.

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
Params.sj(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

% Taxes
Params.tau = 0.15; % Tax rate on labour income
% In addition to payroll tax rate tau, which funds the pension system we will add a progressive 
% income tax which funds government spending.
% The progressive income tax takes the functional form:
% IncomeTax=eta1+eta2*log(Income)*Income; % This functional form is empirically a decent fit for the US tax system
% And is determined by the two parameters
Params.eta1=0.09; % eta1 will be determined in equilibrium to balance gov budget constraint
Params.eta2=0.053;

% Government spending
Params.GdivYtarget = 0.15; % Government spending as a fraction of GDP (this is essentially just used as a target to define a general equilibrium condition)

%% Some initial values/guesses for variables that will be determined in general eqm
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.r=0.1;
Params.AccidentBeq=0.03; % Accidental bequests (this is the lump sum transfer)
Params.G=0.2; % Government expenditure
% Params.eta1=0.09; % already set above

%% Grids
a_grid=10*(linspace(0,1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

% First, z, the AR(1) with age-dependent parameters
[z_grid_J, pi_z_J] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho_z,Params.sigma_epsilon_z,n_z,Params.J);
% z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
% pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j

% Second, e, the iid normal with age-dependent parameters
[e_grid_J, pi_e_J] = discretizeLifeCycleAR1_FellaGallipoliPan(zeros(1,Params.J),Params.sigma_e,vfoptions.n_e,Params.J); % Note: AR(1) with rho=0 is iid normal
% Because e is iid we actually just use
pi_e_J=shiftdim(pi_e_J(1,:,:),1);

% Any (iid) e variable always has to go into vfoptions and simoptions
vfoptions.e_grid=e_grid_J;
vfoptions.pi_e=pi_e_J;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=e_grid_J;
simoptions.pi_e=pi_e_J;


% Grid for labour choice
h_grid=linspace(0,1,n_d)'; % Notice that it is imposing the 0<=h<=1 condition implicitly
% Switch into toolkit notation
d_grid=h_grid;

%% Now, create the return function 
DiscountFactorParamNames={'beta'}; % Just beta, we are using EZ prefs with a warm-glow of bequests, so we had to treat the conditional survival probabilities specially

% Note, 'OLGModel9B_ReturnFn' is just a copy of 'OLGModel6_ReturnFn', except that we have to delete the part about warm-glow of bequests as this has to be treated specially when using EZ. 
% EZ preferences add more risk aversion to this.
ReturnFn=@(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,pension,r,A,delta,alpha,kappa_j,AccidentBeq, eta1,eta2,tau)...
    OLGModel9B_ReturnFn(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,pension,r,A,delta,alpha,kappa_j,AccidentBeq, eta1,eta2,tau);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray'); % Put no households anywhere on grid
jequaloneDist(1,floor((n_z+1)/2),floor((simoptions.n_e+1)/2))=1; % All agents start with zero assets, and the median shock

%% Agents age distribution
% Many OLG models include some kind of population growth, and perhaps
% some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
% Many OLG models include some kind of population growth, and perhaps
% some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one

AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age

%% Test
disp('Test StationaryDist')
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeq','G','eta1'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)
% Note: we need to add z & e to FnsToEvaluate inputs.

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.H = @(h,aprime,a,z,e) h; % Aggregate labour supply
FnsToEvaluate.L = @(h,aprime,a,z,e,kappa_j) kappa_j*exp(z+e)*h;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.K = @(h,aprime,a,z,e) a;% Aggregate  physical capital
FnsToEvaluate.PensionSpending = @(h,aprime,a,z,e,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft = @(h,aprime,a,z,e,sj) aprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.IncomeTaxRevenue = @(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A) OLGModel6_ProgressiveIncomeTaxFn(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.capitalmarket = @(r,K,L,alpha,delta,A) r-alpha*A*(K^(alpha-1))*(L^(1-alpha)); % interest rate equals marginal product of capital net of depreciation
GeneralEqmEqns.pensions = @(PensionSpending,tau,L,r,A,alpha,delta) PensionSpending-tau*(A*(1-alpha)*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*L; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.bequests = @(AccidentalBeqLeft,AccidentBeq,n) AccidentalBeqLeft/(1+n)-AccidentBeq; % Accidental bequests received equal accidental bequests left
GeneralEqmEqns.Gtarget = @(G,GdivYtarget,A,K,L,alpha) G-GdivYtarget*(A*K^(alpha)*(L^(1-alpha))); % G is equal to the target, GdivYtarget*Y
GeneralEqmEqns.govbudget = @(G,IncomeTaxRevenue) G-IncomeTaxRevenue; % Government budget balances (note that pensions are a seperate budget)
% Note: the pensions general eqm condition looks more complicated just because we replaced w with the formula for w in terms of r. It is actually just the same formula as before.

%% Test
% Note: Because we used simoptions we must include this as an input
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid_J,[],simoptions);

%% Solve for the General Equilibrium
heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z_J, d_grid, a_grid, z_grid_J, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm.r;
Params.pension=p_eqm.pension;
Params.AccidentBeq=p_eqm.AccidentBeq;
Params.G=p_eqm.G;
Params.eta1=p_eqm.eta1;

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

figure(1)
subplot(3,1,1); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
title('Life Cycle Profile: Hours Worked')
subplot(3,1,2); plot(1:1:Params.J,AgeConditionalStats.L.Mean)
title('Life Cycle Profile: Labour Supply')
subplot(3,1,3); plot(1:1:Params.J,AgeConditionalStats.K.Mean)
title('Life Cycle Profile: Assets')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel6_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to the FnsToEvaluate
FnsToEvaluate.Consumption=@(h,aprime,a,z,e,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq) OLGModel6_ConsumptionFn(h,aprime,a,z,e,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid_J,[],simoptions);

% GDP
Y=Params.A*(AggVars.K.Mean^Params.alpha)*(AggVars.L.Mean^(1-Params.alpha));

% wage (note that this is calculation is only valid because we have Cobb-Douglas production function and are looking at a stationary general equilibrium)
KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % wage rate (per effective labour unit)

fprintf('Following are some aggregates of the model economy: \n')
fprintf('Output: Y=%8.2f \n',Y)
fprintf('Capital-Output ratio: K/Y=%8.2f \n',AggVars.K.Mean/Y)
fprintf('Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/Y)
fprintf('Average labor productivity: Y/H=%8.2f \n', Y/AggVars.H.Mean)
fprintf('Government-to-Output ratio: G/Y=%8.2f \n', Params.G/Y)
fprintf('Accidental Bequests as fraction of GDP: %8.2f \n',Params.AccidentBeq/Y)
fprintf('Wage: w=%8.2f \n',w)



