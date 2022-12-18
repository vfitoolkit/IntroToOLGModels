%% OLG Model 12: Married Couples
% We model a household as a married couple that has two labor supply
% decisions (one for each spouse), and has two labor productivity processes
% (so a deterministic age-profile, a persistent AR(1) shock, and a
% transitory i.i.d. shock for each spouse).
% 
% We essentially take OLG Model 6, and just change the household problem
% (and so also the functions to evaluate). The state space now has two
% decision variables, one endogenous state (asset holdings), and four
% exogenous states (two markov, z, and two iid, e).

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=[11,11]; % Endogenous labour choice (fraction of time worked)
n_a=301; % Endogenous asset holdings
% Exogenous labor productivity units shocks (next two lines)
n_z=[5,5]; % AR(1) with age-dependent params
vfoptions.n_e=[3,3]; % iid
N_j=Params.J; % Number of periods in finite horizon

vfoptions.lowmemory=1; % The married household problem is larger, so loop over e
simoptions.parallel=3; % Reduce gpu memory requirements (by using sparse cpu matrix instead of gpu matrix; by default it would be parallel=2)

%% Parameters

% Discount rate
Params.beta = 0.96;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 0; % Weight on leisure

Params.A=1; % Aggregate TFP. Not actually used anywhere.
% Production function
Params.alpha = 0.3; % Share of capital
Params.delta = 0.1; % Depreciation rate of capital

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr=46;
% Population growth rate
Params.n=0.02; % percentage rate (expressed as fraction) at which population growths

% Age-dependent labor productivity units (two, one for each spouse)
Params.kappa_j1=[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];
Params.kappa_j2=0.9*[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];
% I am not aware of any age-dependent parameter estimates of correlated
% shocks for married couples. The following are instead loosely based on Wu & Krueger (2021). 
Params.rho_z_married=[0.9,0;0,0.7];
Params.sigmasq_epsilon_z_married=[0.0303, 0.0027; 0.0027, 0.0382]; 
Params.sigma_epsilon_z_married=sqrt(Params.sigmasq_epsilon_z_married);
% The Farmer-Toda method can discretize a VAR(1) with any (postivite semi-definite) variance-covariance matrix.
% iid processes on idiosyncratic labor units which are correlated
Params.sigmasq_epsilon_e_married=[0.1,0.05;0.05,0.1];
Params.sigma_epsilon_e_married=sqrt(Params.sigmasq_epsilon_e_married);

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

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

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
Params.r=0.16;
Params.pension=2; % Initial guess (this will be determined in general eqm)
Params.AccidentBeq=0.07; % Accidental bequests (this is the lump sum transfer) 
Params.G=0.5;  % Government expenditure
Params.eta1=0.51;
% Params.eta1=0.09; % already set above

%% Grids
a_grid=10*(linspace(0,1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

% First, the AR(1) process on z1 and z2 with correlated innovations.
% Notice that this is just a VAR(1) process on z1 and z2 (with zeros on the off-diagonals of the auto-correlation matrix)
% We use the Farmer-Toda method to discretize the VAR(1)
% Note that for VAR(1), the Farmer-Toda method produces a 'joint grid'
[z_grid, pi_z]=discretizeVAR1_FarmerToda([0;0],Params.rho_z_married,Params.sigma_epsilon_z_married,n_z);

% Second, the iid process on e1 and e2 which are correlated
% Notice that this is just a VAR(1) with zero auto-correlation
[e_grid, pi_e]=discretizeVAR1_FarmerToda(zeros(2,1),zeros(2,2),Params.sigma_epsilon_e_married,vfoptions.n_e);
pi_e=pi_e(1,:)';  % Because it is iid, the distribution is just the first row (all rows are identical). We use pi_e as a column vector for VFI Toolkit to handle iid variables.

% To use e variables we need to put them in vfoptions and simoptions
vfoptions.e_grid=e_grid;
vfoptions.pi_e=pi_e;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;
% (Because z_grid does not depend on age we do not need to put it into vfoptions and simoptions)

% Grid for labour choice
h1_grid=linspace(0,1,n_d(1))'; % Notice that it is imposing the 0<=h1<=1 condition implicitly
h2_grid=linspace(0,1,n_d(2))'; % Notice that it is imposing the 0<=h2<=1 condition implicitly
% Switch into toolkit notation
d_grid=[h1_grid; h2_grid];

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

% Notice we use 'OLGModel12_ReturnFn'
ReturnFn=@(h1,h2,aprime,a,z1,z2,e1,e2,sigma,psi,eta,agej,Jr,J,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,warmglow1,warmglow2,AccidentBeq,tau)...
    OLGModel12_ReturnFn(h1,h2,aprime,a,z1,z2,e1,e2,sigma,psi,eta,agej,Jr,J,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,warmglow1,warmglow2,AccidentBeq,tau);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray'); % Put no households anywhere on grid
jequaloneDist(1,floor((n_z(1)+1)/2),floor((n_z(2)+1)/2),floor((simoptions.n_e(1)+1)/2),floor((simoptions.n_e(2)+1)/2))=1; % All agents start with zero assets, and the median shock

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
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeq','G','eta1'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)
% Note: we need to add both of each of the z & e to FnsToEvaluate inputs.

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.H = @(h1,h2,aprime,a,z1,z2,e1,e2) h1+h2; % Aggregate 'hours worked'
FnsToEvaluate.L = @(h1,h2,aprime,a,z1,z2,e1,e2,kappa_j1,kappa_j2) kappa_j1*exp(z1+e1)*h1+kappa_j2*exp(z2+e2)*h2;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.K = @(h1,h2,aprime,a,z1,z2,e1,e2) a;% Aggregate  physical capital
FnsToEvaluate.PensionSpending = @(h1,h2,aprime,a,z1,z2,e1,e2,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft = @(h1,h2,aprime,a,z1,z2,e1,e2,sj) aprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.IncomeTaxRevenue = @(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2) OLGModel12_ProgressiveIncomeTaxFn(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)

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
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);

%% Solve for the General Equilibrium
heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm.r;
Params.pension=p_eqm.pension;
Params.AccidentBeq=p_eqm.AccidentBeq;
Params.G=p_eqm.G;
Params.eta1=p_eqm.eta1;

% Add fraction of time worked for each spouse to FnsToEvaluate
FnsToEvaluate.H1 = @(h1,h2,aprime,a,z1,z2,e1,e2) h1;  % hours worked of first spouse
FnsToEvaluate.H2 = @(h1,h2,aprime,a,z1,z2,e1,e2) h2;  % hours worked of second spouse 

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

% figure(1)
% subplot(4,1,1); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
% title('Life Cycle Profile: Hours Worked')
% subplot(4,1,2); plot(1:1:Params.J,AgeConditionalStats.L.Mean)
% title('Life Cycle Profile: Labour Supply of household')
% subplot(4,1,3); plot(1:1:Params.J,AgeConditionalStats.H1.Mean,1:1:Params.J,AgeConditionalStats.H2.Mean)
% title('Life Cycle Profile: Fraction of time worked, of each spouse')
% legend('First spouse','Second spouse')
% subplot(4,1,4); plot(1:1:Params.J,AgeConditionalStats.K.Mean)
% title('Life Cycle Profile: Assets')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel6_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to the FnsToEvaluate
FnsToEvaluate.Consumption=@(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,tau,AccidentBeq) OLGModel12_ConsumptionFn(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,tau,AccidentBeq);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid,[],simoptions);

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



