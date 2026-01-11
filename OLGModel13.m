%% OLG Model 13: Married Couples, Single Males and Single Females
% We will now essentially combining OLG Models 11 and 12. There are three
% household types: married couples (who make two labor supply decisions,
% and have four exogenous shocks) and single males and single females (one
% labor supply decision and two exogenous shocks; difference between the
% two is in the parameters).
%
% We can set different grid, different ReturnFn, different FnsToEvaluate, etc.,
% for different permanent agent types in the same way as for anything else.
%
% Note that an alterative to having different permanent types use different
% parameters by using 'PTypeNames', that is things like Params.sigma.male=2 and 
% Params.sigma.female=3, you can just use different parameter names for the 
% different permanent types, e.g., Params.sigma_male=2 and
% Params.sigma_female=3. This code mostly uses the PTypeNames for
% differentiating male and female, and uses different parameter names for
% differnetiating married from the other two. This just seemed the easiest
% way to set it up to my mind.

Names_i={'married','male','female'};

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d.married=[11,11]; % Endogenous labour choice (fraction of time worked)
n_d.male=11; % Endogenous labour choice (fraction of time worked)
n_d.female=11; % Endogenous labour choice (fraction of time worked)
n_a=301; % Endogenous asset holdings
% Exogenous labor productivity units shocks (next two lines)
n_z.married=[5,5]; % AR(1) with age-dependent params
n_z.male=5;
n_z.female=5;
vfoptions.n_e.married=[3,3]; % iid
vfoptions.n_e.male=3; % iid
vfoptions.n_e.female=3; % iid
N_j=Params.J; % Number of periods in finite horizon

vfoptions.lowmemory.married=1; % The married household problem is larger, so loop over e
simoptions.parallel.married=3; % Reduce gpu memory requirements (by using sparse cpu matrix instead of gpu matrix; by default it would be parallel=2)
% Because there are a lot of permanent types we will use option that while
% things are calculated on gpu they get stored on cpu. (Otherwise run out of gpu memory)
% Results in minor speed reduction, but managable.
vfoptions.ptypestorecpu=1;
simoptions.ptypestorecpu=1;

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

% Labor productivity process
% Life-cycle AR(1) process z, on (log) labor productivity units
% Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Originals just cover ages 24 to 60, so I create these, and then repeat first and last periods to fill it out
rho_z=0.7596+0.2039*((1:1:37)/10)-0.0535*((1:1:37)/10).^2+0.0028*((1:1:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
sigma_epsilon_z=0.0518-0.0405*((1:1:37)/10)+0.0105*((1:1:37)/10).^2-0.0002*((1:1:37)/10).^3; % Chosen following Karahan & Ozkan (2013) [as used by Fella, Gallipoli & Pan (2019)]
% Note that 37 covers 24 to 60 inclusive
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
rho_z=[rho_z(1)*ones(1,4),rho_z,rho_z(end)*ones(1,4),zeros(1,100-65+1)];
sigma_epsilon_z=[sigma_epsilon_z(1)*ones(1,4),sigma_epsilon_z,sigma_epsilon_z(end)*ones(1,4),sigma_epsilon_z(end)*ones(1,100-65+1)];
% Transitory iid shock
sigma_e=0.0410+0.0221*((24:1:60)/10)-0.0069*((24:1:60)/10).^2+0.0008*((24:1:60)/10).^3;
% Now repeat the first and last values to fill in working age, and put zeros for retirement (where it is anyway irrelevant)
sigma_e=[sigma_e(1)*ones(1,4),sigma_e,sigma_e(end)*ones(1,4),sigma_e(end)*ones(1,100-65+1)];
% Note: These will interact with the endogenous labor so the final labor
% earnings process will not equal that of Karahan & Ozkan (2013)

% Married couple labor productivity process
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
Params.gamma1=0.1; % Fixed-effect
Params.gamma2=-0.1; % Fixed-effect

% Single Female and Single Male labor productivity process
% Age-dependent labor productivity units
Params.kappa_j.female=[linspace(0.5,1.8,Params.Jr-15),linspace(1.8,1,14),zeros(1,Params.J-Params.Jr+1)]; % start the same but peak lower
Params.kappa_j.male=[linspace(0.5,2,Params.Jr-15),linspace(2,1,14),zeros(1,Params.J-Params.Jr+1)];
% Labor productivity shocks
Params.rho_z.female=0.9*rho_z; % Lower autocorrelation (arbitrary, just want to make them different from male)
Params.rho_z.male=rho_z;
Params.sigma_epsilon_z.female=0.9*sigma_epsilon_z; % smaller innovations to the AR(1) (arbitrary, just want to make them different from male)
Params.sigma_epsilon_z.male=sigma_epsilon_z;
Params.sigma_e.male=sigma_e;
Params.sigma_e.female=1.3*sigma_e; % larger iid shocks  (arbitrary, just want to make them different from male)
% Fixed effect
Params.gamma_i.male=0;  % Fixed-effect
Params.gamma_i.female=0;  % Fixed-effect


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
Params.r=0.11;
Params.pension=0.65; % Initial guess (this will be determined in general eqm)
Params.AccidentBeq=0.04; % Accidental bequests (this is the lump sum transfer) 
Params.G=0.2; % Government expenditure
Params.eta1=0.18;

%% Grids
a_grid=10*(linspace(0,1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

% Exogenous shocks for married couple
% First, the AR(1) process on z1 and z2 with correlated innovations.
% Notice that this is just a VAR(1) process on z1 and z2 (with zeros on the off-diagonals of the auto-correlation matrix)
% We use the Farmer-Toda method to discretize the VAR(1)
% Note that for VAR(1), the Farmer-Toda method produces a 'joint grid'
[z_grid, pi_z]=discretizeVAR1_FarmerToda([0;0],Params.rho_z_married,Params.sigma_epsilon_z_married,n_z.married);

% Second, the iid process on e1 and e2 which are correlated
% Notice that this is just a VAR(1) with zero auto-correlation
[e_grid, pi_e]=discretizeVAR1_FarmerToda(zeros(2,1),zeros(2,2),Params.sigma_epsilon_e_married,vfoptions.n_e.married);
pi_e=pi_e(1,:)';  % Because it is iid, the distribution is just the first row (all rows are identical). We use pi_e as a column vector for VFI Toolkit to handle iid variables.

% We will put exogenous shocks for the married couple into the vfoptions and simoptions even though they do not depend on age, just because those
% for the other permanent types do (and we kind of want them all in the same place).
vfoptions.z_grid.married=z_grid;
simoptions.z_grid.married=z_grid;
vfoptions.pi_z.married=pi_z;
simoptions.pi_z.married=pi_z;
simoptions.n_e.married=vfoptions.n_e.married;
vfoptions.e_grid.married=e_grid;
simoptions.e_grid.married=e_grid;
vfoptions.pi_e.married=pi_e;
simoptions.pi_e.married=pi_e;

% Exogenous shocks for single male and single female
% First, z, the AR(1) with age-dependent parameters
[z_grid_J.male, pi_z_J.male] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho_z.male,Params.sigma_epsilon_z.male,n_z.male,Params.J);
[z_grid_J.female, pi_z_J.female] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho_z.female,Params.sigma_epsilon_z.female,n_z.female,Params.J);
% z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
% pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j

% Second, e, the iid normal with age-dependent parameters
[e_grid_J.male, pi_e_J.male] = discretizeLifeCycleAR1_FellaGallipoliPan(zeros(1,Params.J),Params.sigma_e.male,vfoptions.n_e.male,Params.J); % Note: AR(1) with rho=0 is iid normal
[e_grid_J.female, pi_e_J.female] = discretizeLifeCycleAR1_FellaGallipoliPan(zeros(1,Params.J),Params.sigma_e.female,vfoptions.n_e.female,Params.J); % Note: AR(1) with rho=0 is iid normal
% Because e is iid we actually just use
pi_e_J.male=shiftdim(pi_e_J.male(1,:,:),1);
pi_e_J.female=shiftdim(pi_e_J.female(1,:,:),1);

% To use exogenous shocks that depend on age you have to add them to vfoptions and simoptions
vfoptions.z_grid.male=z_grid_J.male;
vfoptions.pi_z.male=pi_z_J.male;
simoptions.z_grid.male=z_grid_J.male;
simoptions.pi_z.male=pi_z_J.male;
vfoptions.z_grid.female=z_grid_J.female;
vfoptions.pi_z.female=pi_z_J.female;
simoptions.z_grid.female=z_grid_J.female;
simoptions.pi_z.female=pi_z_J.female;
% Similarly any (iid) e variable always has to go into vfoptions and simoptions
simoptions.n_e.male=vfoptions.n_e.male;
simoptions.n_e.female=vfoptions.n_e.female;
vfoptions.e_grid.male=e_grid_J.male;
vfoptions.pi_e.male=pi_e_J.male;
simoptions.e_grid.male=e_grid_J.male;
simoptions.pi_e.male=pi_e_J.male;
vfoptions.e_grid.female=e_grid_J.female;
vfoptions.pi_e.female=pi_e_J.female;
simoptions.e_grid.female=e_grid_J.female;
simoptions.pi_e.female=pi_e_J.female;


% Grid for labour choice
h_grid.male=linspace(0,1,n_d.male)'; % Notice that it is imposing the 0<=h<=1 condition implicitly
h_grid.female=linspace(0,1,n_d.female)'; % Notice that it is imposing the 0<=h<=1 condition implicitly
h1_grid=linspace(0,1,n_d.married(1))'; % Notice that it is imposing the 0<=h1<=1 condition implicitly
h2_grid=linspace(0,1,n_d.married(2))'; % Notice that it is imposing the 0<=h2<=1 condition implicitly
% Switch into toolkit notation
d_grid.married=[h1_grid; h2_grid];
% Switch into toolkit notation
d_grid.male=h_grid.male;
d_grid.female=h_grid.female;

% Distribution of the agents across the permanent types (must sum to 1)
Params.ptype_dist=[0.6,0.2,0.2]; % Note: these are same order as Names_i: married, male, female
PTypeDistParamNames={'ptype_dist'};


%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

% Notice we use 'OLGModel12_ReturnFn' for married couple
ReturnFn.married=@(h1,h2,aprime,a,z1,z2,e1,e2,sigma,psi,eta,agej,Jr,J,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,warmglow1,warmglow2,AccidentBeq,tau)...
    OLGModel12_ReturnFn(h1,h2,aprime,a,z1,z2,e1,e2,sigma,psi,eta,agej,Jr,J,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,warmglow1,warmglow2,AccidentBeq,tau);
% And OLGModel10_ReturnFn for the single male and single female
ReturnFn.male=@(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,J,gamma_i,pension,r,A,delta,alpha,kappa_j,warmglow1,warmglow2,AccidentBeq, eta1,eta2,tau)...
    OLGModel10_ReturnFn(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,J,gamma_i,pension,r,A,delta,alpha,kappa_j,warmglow1,warmglow2,AccidentBeq, eta1,eta2,tau);
ReturnFn.female=@(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,J,gamma_i,pension,r,A,delta,alpha,kappa_j,warmglow1,warmglow2,AccidentBeq, eta1,eta2,tau)...
    OLGModel10_ReturnFn(h,aprime,a,z,e,sigma,psi,eta,agej,Jr,J,gamma_i,pension,r,A,delta,alpha,kappa_j,warmglow1,warmglow2,AccidentBeq, eta1,eta2,tau);


%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, vfoptions.z_grid, vfoptions.pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist.married=zeros([n_a,n_z.married,vfoptions.n_e.married],'gpuArray'); % Put no households anywhere on grid
jequaloneDist.married(1,round(n_z.married(1)/2),round(n_z.married(2)/2),round(simoptions.n_e.married(1)/2),round(simoptions.n_e.married(2)/2))=1; % All agents start with zero assets, and the median shock
jequaloneDist.male=zeros([n_a,n_z.male,vfoptions.n_e.male],'gpuArray'); % Put no households anywhere on grid
jequaloneDist.male(1,round(n_z.male/2),round(simoptions.n_e.male/2))=1; % All agents start with zero assets, and the median shock
jequaloneDist.female=zeros([n_a,n_z.female,vfoptions.n_e.female],'gpuArray'); % Put no households anywhere on grid
jequaloneDist.female(1,round(n_z.female/2),round(simoptions.n_e.female/2))=1; % All agents start with zero assets, and the median shock

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
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,simoptions.pi_z,Params,simoptions);

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeq','G','eta1'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)
% Note that because the married and single male/female problems use
% different inputs we need to set up the FnsToEvaluate seperately for each
% permanent type.

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.H.married = @(h1,h2,aprime,a,z1,z2,e1,e2) h1+h2; % Aggregate 'hours worked'
FnsToEvaluate.L.married = @(h1,h2,aprime,a,z1,z2,e1,e2,kappa_j1,kappa_j2) kappa_j1*exp(z1+e1)*h1+kappa_j2*exp(z2+e2)*h2;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.K.married = @(h1,h2,aprime,a,z1,z2,e1,e2) a;% Aggregate  physical capital
FnsToEvaluate.PensionSpending.married = @(h1,h2,aprime,a,z1,z2,e1,e2,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft.married = @(h1,h2,aprime,a,z1,z2,e1,e2,sj) aprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.IncomeTaxRevenue.married = @(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2) OLGModel12_ProgressiveIncomeTaxFn(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)

FnsToEvaluate.H.male = @(h,aprime,a,z,e) h; % Aggregate labour supply
FnsToEvaluate.L.male = @(h,aprime,a,z,e,kappa_j) kappa_j*exp(z+e)*h;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.K.male = @(h,aprime,a,z,e) a;% Aggregate  physical capital
FnsToEvaluate.PensionSpending.male = @(h,aprime,a,z,e,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft.male = @(h,aprime,a,z,e,sj) aprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.IncomeTaxRevenue.male = @(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A,gamma_i,agej,Jr) OLGModel10_ProgressiveIncomeTaxFn(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A,gamma_i,agej,Jr); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)

FnsToEvaluate.H.female = @(h,aprime,a,z,e) h; % Aggregate labour supply
FnsToEvaluate.L.female = @(h,aprime,a,z,e,kappa_j) kappa_j*exp(z+e)*h;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.K.female = @(h,aprime,a,z,e) a;% Aggregate  physical capital
FnsToEvaluate.PensionSpending.female = @(h,aprime,a,z,e,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft.female = @(h,aprime,a,z,e,sj) aprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.IncomeTaxRevenue.female = @(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A,gamma_i,agej,Jr) OLGModel10_ProgressiveIncomeTaxFn(h,aprime,a,z,e,eta1,eta2,kappa_j,r,delta,alpha,A,gamma_i,agej,Jr); % Revenue raised by the progressive income tax (needed own function to avoid log(0) causing problems)


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
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, simoptions.z_grid,simoptions);

%% Solve for the General Equilibrium
heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, [], simoptions.pi_z, d_grid, a_grid, simoptions.z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm.r;
Params.pension=p_eqm.pension;
Params.AccidentBeq=p_eqm.AccidentBeq;
Params.G=p_eqm.G;
Params.eta1=p_eqm.eta1;

% Add fraction of time worked for each spouse to FnsToEvaluate
FnsToEvaluate.H1.married = @(h1,h2,aprime,a,z1,z2,e1,e2) h1;  % hours worked of first spouse
FnsToEvaluate.H2.married = @(h1,h2,aprime,a,z1,z2,e1,e2) h2;  % hours worked of second spouse 
% Notice that it is possible to have FnsToEvaluate that are only for some permanent types but not for others.

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j, Names_i, d_grid, a_grid, vfoptions.z_grid, vfoptions.pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,simoptions.pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,simoptions.z_grid,simoptions);


%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

% figure(1)
% subplot(3,1,1); plot(1:1:Params.J,AgeConditionalStats.H.female.Mean,1:1:Params.J,AgeConditionalStats.H.male.Mean)
% title('Life Cycle Profile: Fraction of time worked, single households')
% legend('single female','single male')
% subplot(3,1,2); plot(1:1:Params.J,AgeConditionalStats.H1.married.Mean,1:1:Params.J,AgeConditionalStats.H2.married.Mean)
% title('Life Cycle Profile: Fraction of time worked, of each spouse of married household')
% legend('First spouse','Second spouse')
% subplot(3,1,3); plot(1:1:Params.J,AgeConditionalStats.K.Mean)
% title('Life Cycle Profile: Assets')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel6_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to the FnsToEvaluate
FnsToEvaluate.Consumption.married=@(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,tau,AccidentBeq) OLGModel12_ConsumptionFn(h1,h2,aprime,a,z1,z2,e1,e2,agej,Jr,pension,r,kappa_j1,kappa_j2,A,alpha,delta,eta1,eta2,tau,AccidentBeq);
FnsToEvaluate.Consumption.male=@(h,aprime,a,z,e,agej,Jr,r,gamma_i,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq) OLGModel10_ConsumptionFn(h,aprime,a,z,e,agej,Jr,r,gamma_i,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq);
FnsToEvaluate.Consumption.female=@(h,aprime,a,z,e,agej,Jr,r,gamma_i,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq) OLGModel10_ConsumptionFn(h,aprime,a,z,e,agej,Jr,r,gamma_i,pension,tau,kappa_j,alpha,delta,A,eta1,eta2,AccidentBeq);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, simoptions.z_grid, simoptions);

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
