%% OLG Model 1: consumption-leisure in general equilibrium OLG
% We start with about the simplest OLG model we can in which there is just the life-cycle with a
% consumption-leisure choice. No assets (savings) and no government.
% The solution will look a bit silly, don't worry, they will improve over the next few models :)
%
% We will solve for the general equilibrium. We then plot some life-cycle
% profiles for the household (based on general equilibrium).
%
% If you have not seen: "Introduction to Life-Cycle Models" I recommend you start there 
% before looking at OLG models
% See: https://www.vfitoolkit.com/updates-blog/2021/an-introduction-to-life-cycle-models/
%
% In this first model households receive pensions which simply 'fall from
% the sky' in the sense that there is no government/tax to fund the
% pensions. This is a bit odd but allows us to focus on the most basic general equilibrium. 
% In the second model we will add a tax to fund pensions.

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
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

% This is the tax which is unused here but will be used to fund the pensions in OLG Model 2
Params.tau=0;
% The amount of pension that retirees receive
Params.pension=0.3;

%% Some initial values/guesses for variables that will be determined in general eqm
Params.w=1;

%% Grids
% While there are no 'a' in this model, VFI Toolkit requires it 
% to figure out what is going on. By making it just a single grid point, 
% and then not using them anywhere, we are essentially solving a model without them.
a_grid=1;
z_grid=[]; % Note: n_z=0 means that z_grid and pi_z will be ignored, but we still need them as inputs to various commands
pi_z=[];

% Grid for labour choice
h_grid=linspace(0,1,n_d)';
% Switch into toolkit notation
d_grid=h_grid;

%% Now, create the return function 
DiscountFactorParamNames={'beta'};

% Note: the following return function actually includes 'assets', but since
% this very basic model has only one possible value for assets this ends up
% irrelevant (VFI Toolkit is not really designed for models this simple with no endogenous state so it looks like it is just overkill when used to solve them). 
ReturnFn=@(h,aprime,a,agej,w,sigma,psi,eta,Jr,pension,tau)...
    OLGModel1_ReturnFn(h,aprime,a,agej,w,sigma,psi,eta,Jr,pension,tau);

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
% This toy model does not have population growth, but there will be
% different fractions of households of each age due to the probability of dying.
% This toy model does not.
Params.mewj=ones(1,Params.J)/Params.J; % This is the Marginal distribution of households over age
% Agents are evenly spread across all ages.

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.

%% Test
disp('Test StationaryDist')
simoptions=struct(); % Just use the defaults
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
toc

%% Everything until this point is just the same as for Life-Cycle Models
% For OLG models we need to add a few things
% 1) GEPriceParamNames: the names of the parameters/prices to be determined in general equilibrium
% 2) FnsToEvaluate: these are same as in Life-Cycle models
% 3) GeneralEqmEqns: the general equilibrium equations, these will be equal to zero in general equilibrium
%               they take any price/parameter as inputs, and also any FnsToEvaluate as inputs
% 4) That is it for set up, now just use HeteroAgentStationaryEqm_Case1_FHorz() command to solve for the general equilibrium
% 5) p_eqm contains the general equilibrium values of parameters/prices determined in general 
%               equilibrium, we can put these in Params and then just calculate anything like value 
%               function, agent stationary distribution, life-cycle profiles, etc., as normal.

%% General eqm variables
GEPriceParamNames={'w'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.L = @(h,aprime,a) h; % Aggregate labour supply
% This version of the model is simple enough that tax revenues and pension expenditures can all just be calculated 'directly'.

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.labormarket = @(w,alpha,L,A) w-(1-alpha)*A*L^(-alpha); % wage equals marginal product of labour
% Note: Inputs to general equilibrium conditions must be either aggregate variables or parameters

%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);

% Advanced: If you did want to test the GeneralEqmEqns to makes sure they are doing what you expect the following commented out lines show how
% % % To be able to test the general equilibrium conditions you need to add the aggregate variables into Params
% % AggVarNames=fieldnames(AggVars);
% % for ii=1:length(AggVarNames)
% %     Params.(AggVarNames{ii})=AggVars.(AggVarNames{ii}).Mean;
% % end
% % GeneralEqmConditionsVec=real(GeneralEqmConditions_Case1_v2(GeneralEqmEqns,Params, 2));

%% Solve for the General Equilibrium
% Use the toolkit to find the equilibrium price index.
% In what follows I use the (default) 'search' approach to calculate the General equilibrium.

heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

%% p_eqm contains the general equilibrium 'prices', to calculate things about the model
% economy in general equilibrium we must put these in Params, and then use them as standard.

% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the general equilibrium
Params.w=p_eqm.w;

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

% We could of course also look at things like life-cycle profiles or
% simulate panel data, but we will wait until later OLG Models to do so.
% E.g., OLG Model 3 plots some life-cycle profiles, OLG Model 4 plots some 
% 'aggregates', like GDP and the capital-output ratio.