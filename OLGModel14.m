%% OLG Model 14: Heterogenous households and heterogeneous firms
% We can essentially just think of households and firms as two different
% permanent types of agents (the household problem is finite horizon, while the 
% firm problem is infinite horizon, but this is fine). Obviously a bunch of
% things related to how we set up the general equilbrium change.

Names_i={'household','firm'};
PTypeDistParamNames={'ptypemass'};
Params.ptypemass=[1,1]; % Mass of households and firms are each equal to one


% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

%% Begin setting up to use VFI Toolkit to solve

% Grid sizes to use for household

% Lets model agents from age 20 to age 100, so 81 periods
Params.agejshifter=19; % Age 20 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % =81, Number of period in life-cycle
n_d.household=51; % Endogenous labour choice (fraction of time worked)
n_a.household=501; % Endogenous share holdings
% Exogenous labor productivity units shocks (next two lines)
n_z.household=15; % AR(1) with age-dependent params
vfoptions.n_e.household=3; % iid
N_j.household=Params.J; % Number of periods in finite horizon

% Grids to use for firm
n_d.firm=101; % Dividend payment
n_a.firm=501; % Capital holdings
n_z.firm=11; % Productivity shock
N_j.firm=Inf; % Infinite horizon

%% Parameters for household

% Discount rate
Params.beta = 1.01; % Changed to get S to increase nearer to 1 given r=0.05 (ran it with beta=0.99, got S=0.3, so increased this; note that it interacts with sj to give the actual discount factor)
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 10; % Weight on leisure

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

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

% Taxes
Params.tau_l = 0.2; % Tax rate on labour income

%% Parameters for firm
% Production
Params.alpha_k=0.311; % diminishing returns to capital input
Params.alpha_l=0.650; % diminishing returns to labor input
Params.delta=0.054; % Deprecitation of physical capital
% Capital adjustment costs
Params.capadjconstant=1.21; % term in the capital adjustment cost
% Tax
Params.tau_corp=0.34; % Tax rate on corporate earnings
Params.phi=0.52; % Fraction of capital adjustment costs that can be deducted from corporate earnings
Params.tau_d=0.2; % Tax rate on dividends
Params.tau_cg=0.2; % Tax rate on capital gains
% Idiosyncatic productivity shocks
Params.rho_z_firm=0.767;
Params.sigma_z_e_firm=0.211;

% Set the firm discount factor below (as it is determined in general eqm)
% Params.firmbeta=1/(1+Params.r/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax

%% Also, we want to target a capital-output ratio
% This is relevant to the general equilibrium conditions
Params.TargetKdivL=2.03;

%% Remaining parameters

% Some initial values/guesses for variables that will be determined in general eqm
Params.pension=0.4; % Initial guess (this will be determined in general eqm)
Params.r=0.05;
Params.w=1;
Params.AccidentBeq=0.03; % Accidental bequests (this is the lump sum transfer)
Params.G=0.1; % Government expenditure
Params.firmbeta=1/(1+Params.r/(1-Params.tau_cg)); % 1/(1+r) but returns net of capital gains tax
Params.D=0.2; % Dividends
Params.P0=1;
Params.Lhscale=0.3; % Scaling the household labor supply

%% Grids for household

% Grid for labour choice
d_grid.household=linspace(0,1,n_d.household)'; % Notice that it is imposing the 0<=h<=1 condition implicitly

% Grid for share holdings
a_grid.household=10*(linspace(0,1,n_a.household).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

% First, z, the AR(1) with age-dependent parameters
[z_grid_J, pi_z_J] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho_z,Params.sigma_epsilon_z,n_z.household,Params.J);
% z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
% pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j

% Second, e, the iid normal with age-dependent parameters
[e_grid_J, pi_e_J] = discretizeLifeCycleAR1_FellaGallipoliPan(zeros(1,Params.J),Params.sigma_e,vfoptions.n_e.household,Params.J); % Note: AR(1) with rho=0 is iid normal
% Because e is iid we actually just use
pi_e_J=shiftdim(pi_e_J(1,:,:),1);

% z_grid and pi_z for household
z_grid.household=z_grid_J;
pi_z.household=pi_z_J;

% Any (iid) e variable always has to go into vfoptions and simoptions
vfoptions.e_grid.household=e_grid_J;
vfoptions.pi_e.household=pi_e_J;
simoptions.n_e.household=vfoptions.n_e.household;
simoptions.e_grid.household=e_grid_J;
simoptions.pi_e.household=pi_e_J;


%% Grids for firm
d_grid.firm=linspace(0,1,n_d.firm)'; % Notice that it is imposing the d>=0 condition implicitly
a_grid.firm=10*(linspace(0,1,n_a.firm).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.

[z_grid.firm,pi_z.firm] = discretizeAR1_FarmerToda(0,Params.rho_z_firm,Params.sigma_z_e_firm,n_z.firm);
z_grid.firm=exp(z_grid.firm);


%% Now, create the return function

% For households
DiscountFactorParamNames.household={'beta','sj'};
% Notice we use 'OLGModel14_HouseholdReturnFn'
ReturnFn.household=@(h,sprime,s,z,e,sigma,psi,eta,agej,Jr,J,pension,r,w,P0,D,kappa_j,warmglow1,warmglow2,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale)...
    OLGModel14_HouseholdReturnFn(h,sprime,s,z,e,sigma,psi,eta,agej,Jr,J,pension,r,w,P0,D,kappa_j,warmglow1,warmglow2,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale);

% For firms
DiscountFactorParamNames.firm={'firmbeta'};
% Notice we use 'OLGModel14_FirmReturnFn'
ReturnFn.firm=@(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg)...
    OLGModel14_FirmReturnFn(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
disp('Test ValueFnIter')
tic;
% Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
toc

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist.household=zeros([n_a.household,n_z.household,vfoptions.n_e.household],'gpuArray'); % Put no households anywhere on grid
jequaloneDist.household(1,floor((n_z.household+1)/2),floor((simoptions.n_e.household+1)/2))=1; % All agents start with zero assets, and the median shock

% Note that because the firms are infinite horizon they do not have an age=1 distribution

%% Agents age distribution
% Many OLG models include some kind of population growth, and perhaps some other things that create a weighting of different ages that needs to
% be used to calculate the stationary distribution and aggregate variable.
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one

AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age

%% Test
disp('Test StationaryDist')
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% General eqm variables
GEPriceParamNames={'r','pension','AccidentBeq','G','w','firmbeta','D','P0','Lhscale'}; 
% We don't need P
% We can get P from the equation that defines r as the return to the mutual fund
% 1+r = (P0 +(1-tau_d)D - tau_cg(P0-P))/Plag
% We are looking at stationary general eqm, so
% Plag=P;
% And thus we have
% P=((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg);

%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)
% Note: we need to add z & e to FnsToEvaluate inputs for households,
% whereas firm only has z (it is just coincidence/lazy that I call them
% both z).

% Stationary Distribution Aggregates from households (important that ordering of Names and Functions is the same)
FnsToEvaluate.L_h.household = @(h,sprime,s,z,e,kappa_j,Lhscale) kappa_j*exp(z+e)*Lhscale*h;  % Aggregate labour supply in efficiency units 
FnsToEvaluate.S.household = @(h,sprime,s,z,e) s; % Aggregate share holdings
FnsToEvaluate.PensionSpending.household = @(h,sprime,s,z,e,pension,agej,Jr) (agej>=Jr)*pension; % Total spending on pensions
FnsToEvaluate.PayrollTaxRevenue.household = @(h,sprime,s,z,e,agej,Jr,tau_l,w,kappa_j,Lhscale) (agej<Jr)*tau_l*w*kappa_j*exp(z+e)*Lhscale*h; % Total spending on pensions
FnsToEvaluate.AccidentalBeqLeft.household = @(h,sprime,s,z,e,sj) sprime*(1-sj); % Accidental bequests left by people who die
FnsToEvaluate.CapitalGainsTaxRevenue.household = @(h,sprime,s,z,e,tau_cg,P0,D,tau_d,r) tau_cg*(P0-(((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg)))*s; % tau_cg*(P0-Plag)*s, but substitute P=Plag, and then substitute for P
% From firms
FnsToEvaluate.Output.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) z*(k^alpha_k)*((w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)))^alpha_l; % Production function z*(k^alpha_k)*(l^alpha_l) (substituting for l)
FnsToEvaluate.L_f.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) (w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % (effective units of) labor demanded by firm
FnsToEvaluate.K.firm = @(d,kprime,k,z,w,alpha_k,alpha_l) k; % physical capital
FnsToEvaluate.DividendPaid.firm = @(d,kprime,k,z,w) d; % dividend paid by firm
FnsToEvaluate.Sissued.firm = @(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) OLGModel14_FirmShareIssuance(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % Share issuance
FnsToEvaluate.CorpTaxRevenue.firm = @(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi) OLGModel14_FirmCorporateTaxRevenue(d,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi); % revenue from the corporate profits tax

% General Equilibrium conditions (these should evaluate to zero in general equilbrium)
GeneralEqmEqns.sharemarket = @(S) S-1; % mass of all shares equals one
GeneralEqmEqns.labormarket = @(L_h,L_f) L_h-L_f; % labor supply of households equals labor demand of firms
GeneralEqmEqns.pensions = @(PensionSpending,PayrollTaxRevenue) PensionSpending-PayrollTaxRevenue; % Retirement benefits equal Payroll tax revenue: pension*fractionretired-tau*w*H
GeneralEqmEqns.bequests = @(AccidentalBeqLeft,AccidentBeq,n) AccidentalBeqLeft/(1+n)-AccidentBeq; % Accidental bequests received equal accidental bequests left
GeneralEqmEqns.govbudget = @(G,tau_d,D,CapitalGainsTaxRevenue,CorpTaxRevenue) G-tau_d*D-CapitalGainsTaxRevenue-CorpTaxRevenue; % G is equal to the target, GdivYtarget*Y
GeneralEqmEqns.firmdiscounting = @(firmbeta,r,tau_cg) firmbeta-1/(1+Params.r/(1-Params.tau_cg)); % Firms discount rate is related to market return rate
GeneralEqmEqns.dividends = @(D,DividendPaid) D-DividendPaid; % That the dividend households receive equals that which firms give
GeneralEqmEqns.ShareIssuance = @(Sissued,P0,D,tau_cg,tau_d,r) P0-((((1-tau_cg)*P0 + (1-tau_d)*D)/(1+r-tau_cg))-Sissued); % P0=P-S, but substitute for P (see derivation inside the return fn)
GeneralEqmEqns.CapitalOutputRatio =@(K,L_f,TargetKdivL) K/L_f-TargetKdivL;

%% Test
% Note: Because we used simoptions we must include this as an input
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j,Names_i, d_grid, a_grid, z_grid,simoptions);

% Next few lines were used to try a few parameter values so as to get a
% decent initial guess before actually solving the general equilbrium
fprintf('Check: L_h, L_f, K \n')
[AggVars.L_h.Mean,AggVars.L_f.Mean,AggVars.K.Mean]
fprintf('Check: K/L_f (should be about 2.03 \n')
AggVars.K.Mean/AggVars.L_f.Mean
fprintf('Check: S \n')
AggVars.S.Mean
fprintf('Check: ShareIssuance GE condition \n')
Params.P0-((((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg))-AggVars.S.Mean)


%% Solve for the General Equilibrium
% heteroagentoptions.fminalgo=4 % CMA-ES algorithm 

heteroagentoptions.verbose=1;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% p_eqm contains the general equilibrium parameter values
% Put this into Params so we can calculate things about the initial equilibrium
Params.r=p_eqm.r;
Params.pension=p_eqm.pension;
Params.AccidentBeq=p_eqm.AccidentBeq;
Params.G=p_eqm.G;
Params.w=p_eqm.w;
Params.firmbeta=p_eqm.firmbeta;
Params.D=p_eqm.D;
Params.P0=p_eqm.P0;
Params.Lhscale=p_eqm.Lhscale;

% Calculate a few things related to the general equilibrium.
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j, Names_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
% Can just use the same FnsToEvaluate as before.
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Plot the life cycle profiles of capital and labour for the inital and final eqm.

figure(1)
subplot(2,1,1); plot(1:1:Params.J,AgeConditionalStats.L_h.Mean)
title('Life Cycle Profile: Effective Labour Supply')
subplot(2,1,2); plot(1:1:Params.J,AgeConditionalStats.S.Mean)
title('Life Cycle Profile: Share holdings')
% saveas(figure_c,'./SavedOutput/Graphs/OLGModel6_LifeCycleProfiles','pdf')

%% Calculate some aggregates and print findings about them

% Add consumption to the FnsToEvaluate
FnsToEvaluate.Consumption.household=@(h,sprime,s,z,e,agej,Jr,pension,r,w,P0,D,kappa_j,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale) OLGModel14_HouseholdConsumptionFn(h,sprime,s,z,e,agej,Jr,pension,r,w,P0,D,kappa_j,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale);
FnsToEvaluate.Income.household=@(h,sprime,s,z,e,agej,Jr,pension,r,w,P0,D,kappa_j,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale) OLGModel14_HouseholdIncomeFn(h,sprime,s,z,e,agej,Jr,pension,r,w,P0,D,kappa_j,AccidentBeq,tau_l,tau_d,tau_cg,Lhscale);

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid,simoptions);

Y=AggVars.Output.Mean;

P=((1-Params.tau_cg)*Params.P0 + (1-Params.tau_d)*Params.D)/(1+Params.r-Params.tau_cg);

% Calculate the aggregate TFP as output/((capital^alpha_k)*(labor^alpha_l))
AggregateTFP=Y/((AggVars.K.Mean^Params.alpha_k)*(AggVars.L_f.Mean^Params.alpha_l));

% Total value of firms
temp=V.firm.*StationaryDist.firm;
temp(StationaryDist.firm==0)=0; % Get rid of points that have V=-inf but zero mass which would give nan
TotalValueOfFirms=sum(sum(temp));

fprintf('Following are some aggregates of the model economy: \n')
fprintf('Output: Y=%8.2f \n',AggVars.Output.Mean)
fprintf('Aggregate TFP: Y=%8.2f \n',AggregateTFP)
fprintf('Capital-Output ratio (firm side): K/Y=%8.2f \n',AggVars.K.Mean/Y)
fprintf('Total asset value (HH side): P*S=%8.2f \n',P*AggVars.S.Mean)
fprintf('Total firm value (firm side): Value of firm=%8.2f \n',TotalValueOfFirms)
fprintf('Consumption-Output ratio: C/Y=%8.2f \n',AggVars.Consumption.Mean/Y)
fprintf('Government-to-Output ratio: G/Y=%8.2f \n', Params.G/Y)
fprintf('Wage: w=%8.2f \n',Params.w)










