%% OLG Model 3: Dependence on Age
% We will just extend the model in three simple ways. We will let wages have
% a 'life cycle profile'. We will include 'conditional survival
% probabilities', so as people get older their probability of surviving one
% more year (until next period) falls. This will give our model a
% population an age-demographic that looks more realistic.
% To further add realism to the age-demographics we will also add constant
% population growth rate n.
% We will look at two 'new' outputs: life-cycle profile, and demographic pyramid
% 
% Note that this is not really changing anything much about general
% equilibrium. It is more about adding a bit more realistic life-cycle.
% 
% Just to show how easy it is I am also going to change so that instead of
% fixing pensions and determining tax in general equilibrium (to balance
% the budget of the pension system) we will instead fix the tax and
% determine the pensions in general equilibrium.
% Note that this essentially just involves changing GEPriceParamNames
% 
% Comment: Normally when adding conditional survival probabilities you need
% to take care about people leaving behind assets when they die, but this
% model currently has no assets so it is not an issue yet.

%% Changes from OLG Model 2
% Add Params.kappa_j, the labor productivity units conditional on age
% Add Params.sj, the survival probabilities
% Adjust Params.mewj, based on Params.sj
% Draw Life-Cycle Profiles
% Draw demographic-pyramid
% Modify FnsToEvaluate to give both labor supply and hours worked

% Note that both sj and kappa_j are vectors that depend on age j. VFI
% Toolkit detects this automatically and deals with them appropriately in
% all the commands. (It detects that length() of the parameter shows a
% vector of length N_j; which is set to be the number of periods/ages).

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Initial age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=101; % Endogenous labour choice (fraction of time worked)
n_a=1;
n_z=0; % This is how the VFI Toolkit thinks about deterministic models
N_j=Params.J; % Number of periods in finite horizon

figure_c=0; % I like to use a counter for the figures. Makes it easier to keep track of them when editing.

%% Parameters

% Discount rate
Params.beta = 0.99;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 10; % Weight on leisure

Params.A=1; % Aggregate TFP. Not actually used anywhere.
% Production function
Params.alpha = 0.3; % Share of capital
Params.delta = 0.1; % Depreciation rate of capital

% Demographics
Params.Jr=67-Params.agejshifter; % Retirement age is 67 (remember j=1 is age 20) (You work when younger, not working from Jr on)
Params.agej=(1:1:Params.J)'; % Current 'j' age, so can check when you retire.
% Population growth rate
Params.n=0.02; % percentage rate (expressed as fraction) at which population growths

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

% Labor efficiency units depend on age
Params.kappa_j=[linspace(1,3,50-Params.agejshifter), linspace(3,2,(Params.Jr-1)-(50-Params.agejshifter)),zeros(1,Params.J-Params.Jr+1)];
    % These are not done seriously, really they should be set to something like 'average hourly wage conditional on age' in the data.
    % I have made them increase until age 50 (j=31), then decrease, and then be zero from retirement at age 67.
    % For parameters that depend on age, we just make them a vector with a
    % length the same as the number of periods, VFI Toolkit then handles
    % them automatically.

% Taxes
Params.tau = 0.15; % Tax rate on labour income

%% Some initial values/guesses for variables that will be determined in general eqm
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.w=1;

%% Grids
% While there are no 'a' for 'z' in this model, VFI Toolkit requires them 
% to figure out what is going on. By making them just a single grid point, 
% and then not using them anywhere, we are essentially solving a model without them.
a_grid=1;
z_grid=1;  % Note: n_z=0 means that z_grid and pi_z will be ignored, but we still need them as inputs to various commands
pi_z=1;

% Grid for labour choice
h_grid=linspace(0,1,n_d)';
% Switch into toolkit notation
d_grid=h_grid;

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

% We need to change the return function so that it includes kappa_j (compared to FiscalOLG1)
ReturnFn=@(h,aprime,a,agej,w,sigma,psi,eta,Jr,pension,tau,kappa_j)...
    OLGModel3_ReturnFn(h,aprime,a,agej,w,sigma,psi,eta,Jr,pension,tau,kappa_j);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
vfoptions=struct(); % Just using the defaults.
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
jequaloneDist=1; % n_a by n_z, but in current setup this is just 1-by-1 in anycase so 'everyone' is going to have to be born here.

%% Agents age distribution
% Many OLG models include some kind of population growth, and perhaps
% some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

%% Test
disp('Test StationaryDist')
simoptions=struct(); % Just use the defaults
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
toc

%% Plot the demographic pyramid
% Note that because aging is completely exogenous it is not affected by the
% general equilibrium (initial or final), nor does it change over the
% transition. (We will allow this in a more advanced model later on.)
figure_c=figure_c+1;
figure(figure_c)
hold on
pyramidright = barh((1:1:Params.J)+Params.agejshifter,Params.mewj/2,'hist'); % I need the divided by two to get the symmetry (half on each 'side')
pyramidleft = barh((1:1:Params.J)+Params.agejshifter,-Params.mewj/2,'hist'); % 'minus', so it goes to other 'side'
set(pyramidright,'FaceColor','b')
set(pyramidleft,'FaceColor','b')
hold off
ylabel('Age')
xlabel('Fraction of population')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel3_DemographicPyramid','pdf')

%% General eqm variables
GEPriceParamNames={'w','pension'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.L = @(h, aprime,a, kappa_j) kappa_j*h; % Aggregate labour supply
FnsToEvaluate.H = @(h, aprime,a) h; % Aggregate hours worked
FnsToEvaluate.PensionSpending = @(h,aprime,a,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
% This version of the model is simple enough that tax revenues can just be calculated 'directly'.
% Because of the surivival probabilities it is no longer as trivial to
% calculate pension spending. We could essentially pre-calculate it (it only depends on
% mewj and on pension), but because of how the toolkit works it is just as simple to just get the code to do it for us.

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.labormarket = @(w,alpha,L,A) w-(1-alpha)*A*L^(-alpha); % wage equals marginal product of labour
GeneralEqmEqns.pensions = @(PensionSpending,tau,w,L) PensionSpending-tau*w*L; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H

%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);

%% Solve for the General Equilibrium

heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the general equilibrium
Params.w=p_eqm.w;
Params.pension=p_eqm.pension;

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid);

%% Plot the life cycle profiles of hours worked and effective labour supply in the general eqm.
figure_c=figure_c+1;
figure(figure_c)
subplot(2,1,1); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
title('Life Cycle Profile: Hours Worked')
subplot(2,1,2); plot(1:1:Params.J,AgeConditionalStats.L.Mean)
title('Life Cycle Profile: Labour Supply')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel3_LifeCycleProfiles','pdf')
