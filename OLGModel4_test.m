%% OLG Model 4: Assets and Aggregate Physical Capital
% We add assets/savings. This involves adding household assets as an endogenous state.
% To fit empirical fact that households often die with substantial asset
% holdings we will include a 'warm glow of bequests' (utility from having
% assets at death)
% Still no government, just the pension system.
%
% When we add assets we now have a representative firm with Cobb-Douglas
% production function, Y=A*(K^alpha)*(L^(1-alpha));
% Together with competitive labor and capital markets this means there is a
% tight relationship (an isomorphic function) between r and w. We take
% advantage of this to solve the model in terms of r, and get a formula for
% w in terms of r as follows:
% Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
% KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
% From K/Y, substituting the production function gives
% KdivY=(KdivL^(1-alpha))/A; % Don't use this, just including it for completeness. 
% Notice that discount factor beta plays an indirect role in determining K/Y via determining r
% We know w=MPL (MPL is marginal product of labour)
% w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)
% [You could solve model in terms of w, and use inverse of this formula to
% get r. We simply follow tradition and solve in terms of r. It is easier to 
% get analytical initial guesses for r than for w.]
%
% Because people die unexpectedly (due to sj, the conditional surival
% probabilities) some people will still have assets when they die. The warm
% glow bequest is another reason why people will still have assets when
% they die. We don't want these assets to simply disappear. We will
% therefore redistribute them as a lump-sum to everyone, which we will call
% accidental bequests. This involves adding them to the return function,
% adding them as a function to evaluate, and adding a general equilibrium
% condition that the accidental bequest people leave when they die are
% equal to the received lump-sum transfers.
%

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 20 to age 100, so 81 periods

Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle

% Grid sizes to use
n_d=101; % Endogenous labour choice (fraction of time worked)
n_a=301;
n_z=1; % This is how the VFI Toolkit thinks about deterministic models
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

% Warm-glow of bequest
Params.warmglowparam1=1;
Params.warmglowparam2=2;

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
Params.sj(end)=0; % We did not need this in the previous model, but now it is important to help keep track of accidental bequests (otherwise we would have to use a more complicated formula)

% Labor efficiency units depend on age
Params.kappa_j=[linspace(1,3,50-Params.agejshifter), linspace(3,2,(Params.Jr-1)-(50-Params.agejshifter)),zeros(1,Params.J-Params.Jr+1)];
    % These are not done seriously, really they should be set to something like 'average hourly wage conditional on age' in the data.
    % I have made them increase until age 50 (j=31), then decrease, and then be zero from retirement at age 67.
    % For parameters that depend on age, we just make them a vector with a
    % length the same as the number of periods, VFI Toolkit then handles
    % them automatically.

% Taxes
Params.tau = 0.15; % Tax rate on labour income
Params.tau_final=0.1; % We will look at transition path resulting from reduction in tax rate to this level
Params.tau=Params.tau;


%% Some initial values/guesses for variables that will be determined in general eqm
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.r=0.03;
Params.AccidentBeq=0.03; % Accidental bequests (this is the lump sum transfer)

%% Grids
Params.amax=5; % There is no easy way to figure out what amax should be. 
 % It needs to be large enough that no agent would ever choose to go there 
 % and stay there (you can use policy function to look at this). But we
 % don't want it larger than necessary, as this would 'waste' grid points
 % and so be less accurate. Because the current OLG has no 'idiosycratic
 % shocks' we don't need amax to be very high (relative to income). Later,
 % we will add idiosyncratic shocks and need to set a higher amax.
a_grid=Params.amax*(linspace(0,1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
% Note: Implicitly, we are imposing borrowing constraint a'>=0
z_grid=1;
pi_z=1;

% Grid for labour choice
h_grid=linspace(0,1,n_d)';
% Switch into toolkit notation
d_grid=h_grid;

%%
DiscountFactorParamNames={'beta','sj'};

% We need to change the return function so that it uses r and includes assets
ReturnFn=@(h,aprime,a,z,agej,r,A,delta,alpha,sigma,psi,eta,Jr,pension,tau,kappa_j,J,warmglowparam1,warmglowparam2,AccidentBeq) OLGModel4_ReturnFn(h,aprime,a,z,agej,r,A,delta,alpha,sigma,psi,eta,Jr,pension,tau,kappa_j,J,warmglowparam1,warmglowparam2, AccidentBeq)

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
vfoptions=struct(); % Just using the defaults.
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
jequaloneDist=zeros(n_a,n_z); % n_a by n_z.
jequaloneDist(1,1)=1; % Everyone is born with zero assets

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

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeq'};

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Stationary Distribution Aggregates (important that ordering of Names and Functions is the same)
FnsToEvaluate.H = @(h,aprime,a,z) h; % Aggregate labour supply
FnsToEvaluate.L = @(h,aprime,a,z,kappa_j) kappa_j*h; % Aggregate labour supply (in efficiency units)
FnsToEvaluate.K = @(h,aprime,a,z) a; % Aggregate physical capital
FnsToEvaluate.PensionSpending = @(h,aprime,a,z,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft = @(h,aprime,a,z,sj) aprime*(1-sj); % Accidental bequests left by people who die

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.capitalmarket = @(r,K,L,alpha,delta,A) r-alpha*A*(K^(alpha-1))*(L^(1-alpha)); % interest rate equals marginal product of capital net of depreciation
GeneralEqmEqns.pensions = @(PensionSpending,tau,L,r,A,alpha,delta) PensionSpending-tau*(A*(1-alpha)*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*L; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.bequests = @(AccidentalBeqLeft,AccidentBeq,n) AccidentalBeqLeft/(1+n)-AccidentBeq; % Accidental bequests received equal accidental bequests left
% Note: the pensions general eqm condition looks more complicated just because we replaced w with the formula for w in terms of r. It is actually just the same formula as before.


%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);

%% Solve for the General Equilibrium
heteroagentoptions.fminalgo=5 % 5: 
heteroagentoptions.fminalgo5.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'capitalmarket','r',0,0.2,... % capitalmarket GE condition will be positive if r is too big, so subtract
    'pensions','pension',0,0.2,... % pensions GE condition will be positive if pension is too big, so subtract
    'bequests','AccidentBeq',1,0.2,... % bequests GE condition will be positive if AccidentBeq is too small, so add
    };
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it when add=0
% A small 'factor' will make the convergence to solution take longer, but
% too large a value will make it unstable (fail to converge). Technically
% this is the damping factor in a shooting algorithm.

heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,n_d, n_a, n_z, N_j, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm.r;
Params.pension=p_eqm.pension;
Params.AccidentBeq=p_eqm.AccidentBeq;

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid);

%% Plot the life cycle profiles of capital and labour for the general eqm.

figure_c=figure_c+1;
figure(figure_c)
subplot(3,1,1); plot(1:1:Params.J,AgeConditionalStats.H.Mean)
title('Life Cycle Profile: Hours Worked')
subplot(3,1,2); plot(1:1:Params.J,AgeConditionalStats.L.Mean)
title('Life Cycle Profile: Labour Supply')
subplot(3,1,3); plot(1:1:Params.J,AgeConditionalStats.K.Mean)
title('Life Cycle Profile: Assets')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel4_LifeCycleProfiles','pdf')

%% Often in OLG models we are interested in things like the aggregate capital-output (a.k.a, wealth-income)
% ratio. Let's calculate some aggregates and ratios of these.

% Add consumption to the FnsToEvaluate
FnsToEvaluate.Consumption=@(h,aprime,a,z,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A) OLGModel4_ConsumptionFn(h,aprime,a,z,agej,Jr,r,pension,tau,kappa_j,alpha,delta,A);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z,N_j, d_grid, a_grid, z_grid);

% GDP
Y=Params.A*(AggVars.K.Mean^Params.alpha)*(AggVars.L.Mean^(1-Params.alpha));

fprintf('Following are some aggregates of the model economy: \n')
fprintf('Output: Y=%8.2f (this number is fairly meaningless in and of itself) \n',Y)
fprintf('Capital-Output ratio: K/Y=%8.2f \n',AggVars.K.Mean/Y)
fprintf('Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/Y)
fprintf('Average labor productivity: Y/H=%8.2f \n', Y/AggVars.H.Mean)
fprintf('Accidental Bequests as fraction of GDP: %8.2f \n',Params.AccidentBeq/Y)











