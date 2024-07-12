function sol=...
    optimize_storage_control_for_profit(P_demand, P_re, P_cap_gas, P_cap_storage, cost_to_store, LMP, eta_storage)

% Determine P_curtail, P_gas, P_storage_in, P_storage_out, E_storage, storage_losses, E_cap.
% This data is stored in the outputted data struct, sol.

% Follow MathWorks optimization examples here:
% https://www.mathworks.com/help/optim/examples.html?s_cid=doc_ftr
% Especially this example:
% https://www.mathworks.com/help/optim/ug/optimal-dispatch-problem-based.html

% Limit rate of energy flow into/out of storage to be P_farm

% Uncomment to use fixed inputs:
% load P_demand, P_re, P_cap_gas, P_farm, E_cap, cost_to_store, LMP, eta_storage
% load('input_data_test_optimization.mat');

% NOTE OCT. 7, 2020: In code below, there are some if-statements for if
% P_cap_gas > 0. These have been disabled since that seems to open up many
% new questions and challending sizing constraints. For our current paper,
% we will not have an special conditions on P_cap_gas, even when
% RE_penetration = 100%. That is, even when RE_penetration= 100%, we will
% simualte a scenario where there is a gas plant available with a capacity
% equal to the peak demand.

% NOTE JAN 20, 2021: removing E_stored dependency from cost function.

% NOTE FEB 3, 2021: Updated Van den Bergh costs to be in 2015 Euros converted to 2020 USD
euro_to_usd_2015 = 1.11;
inflation_2015_to_2020 = 1.0342; 

plotOn= 0; % 1- plot resuts. 0- don't plot results.

%% Initialize optimization 

P_demand_year= P_demand';
P_re_year= P_re';
LMP_year = LMP';

% Initialize outputs
sol.P_storage_in=    [];
sol.P_storage_out=   [];
sol.E_storage=       [];
sol.P_curtail=       [];
sol.P_gas=           [];
sol.gas_on=          [];
sol.gas_turns_on=    [];
sol.P_gas_ramp_up=   [];
sol.P_gas_ramp_down= [];
sol.exitflag=        1; % Initialize as no issues

 % break the annual optimization into nPeriods for independent optimization
 % Use final values from previous step for initial conditions 
 % This probably has negligible effect, since optimizaiton effects are probably
 % more in the day-to-day range.
 % Optimizing the entire year at once stalls the optimization solver
nPeriods= 12;

length_separate_periods= floor( length(P_demand) / nPeriods);
inds_start= 1: length_separate_periods: length(P_demand);
inds_start(end)= [];

inds_end=  length_separate_periods : length_separate_periods: length(P_demand);
inds_end(end)= length(P_demand);

for i= 1:nPeriods
           
    P_demand= P_demand_year( inds_start(i) : inds_end(i) );
    P_re= P_re_year( inds_start(i) : inds_end(i) )';
    LMP = LMP_year( inds_start(i) : inds_end(i) );

    num_hours= length(P_demand);

    % Gas plant operation constraints From Kumar paper
    % Cost every time start gas plant from off:
    % cost_to_start_Gas= .24 * P_cap_gas; 
    % This is not used because adding an integer flag for turning on gas plant
    % seemed to stall the optimization.
    % It is still a reasonable optimization becuase we still add the
    % restriction that when the Gas plant is on, it must produce more than a
    % P_min_gas power amount. This is a realistic constraint. The
    % cost_to_start_Gas constraint might make the optimization more optimal,
    % but I would estimate it is not a huge impact, and maybe relatively
    % constrainted already for power balance.

    P_min_gas= 0.2 * P_cap_gas;


    linprob = optimproblem('ObjectiveSense','maximize');

    % Initialize variables
    P_storage_in= optimvar(   'P_storage_in',    num_hours, 'LowerBound', 0, 'UpperBound', P_cap_storage);
    P_storage_out= optimvar(  'P_storage_out',   num_hours, 'LowerBound', 0, 'UpperBound', P_cap_storage);
    E_storage= optimvar(      'E_storage',       num_hours, 'LowerBound',0);%,  'UpperBound',E_cap); 
    P_curtail= optimvar(      'P_curtail',       num_hours, 'LowerBound',0); 
%     if P_cap_gas > 0
    gas_on= optimvar(         'gas_on',          num_hours, 'Type','integer','LowerBound', 0, 'UpperBound',1);
    gas_turns_on= optimvar(   'gas_turns_on',    num_hours, 'Type','integer','LowerBound', 0, 'UpperBound',1);
    P_gas= optimvar(          'P_gas',           num_hours, 'LowerBound',0,  'UpperBound', P_cap_gas);
    P_gas_ramp_up= optimvar(  'P_gas_ramp_up',   num_hours, 'LowerBound', 0);
    P_gas_ramp_down= optimvar('P_gas_ramp_down', num_hours, 'LowerBound', 0);
%     end

    % P_curtail is Power lost when need to shut down farm due to oversupply
    % Van den Bergh applied a cost around LMP 
    % Van den Bergh uses the costs below for gas plant ramping costs
    % cost_ramp= 1; % $ per change in MW 
    % cost_start= 50*P_min_gas; % from Van den Bergh PDF page 20
    % cost_curtail= 10000; % cost per MWh. From Van den Bergh.

    % MNH updated 2021-02-03 - cost values from Van den Bergh in 2015 euros. Adjusted so they match reference and 
    % converted to 2015 USD and adjusted for inflation to 2020 
    % Van den Bergh: table on page 19 in the Van den Bergh paper shows that direct start and indirect start costs sum to 42.4 euros per MW change.
    cost_ramp= 0.8*euro_to_usd_2015*inflation_2015_to_2020; % $ per change in MW
    cost_start = 42.4*euro_to_usd_2015*inflation_2015_to_2020*P_min_gas; % every time the plant is "started up", it needs to go from 0 MW to "p_min" MW 
    cost_curtail= 10000*euro_to_usd_2015*inflation_2015_to_2020; % cost per MWh. From Van den Bergh, 2015 Euros to 2020 USD


    % A numerical value so optimization curtails as last resort

    % NOTE: P_storage_in is power that leaves grid when entering storage
    %       (before input storage inefficiency)
    %       P_storage_out is power that re-enters grid after leaving storage
    %       (after storage output inefficiency)   
    % storage power losses=  (1-eta_storage)*P_storage_in               % P_loss enter storage
    %                      + (1-eta_storage)/eta_storage* P_storage_out % P_loss leave storage
    % total_storage_Energy_losses= sum(storage power losses) 
    %                            = (1-eta_storage^2)*P_storage_in

    %% Define constraints

    if isempty( sol.E_storage )
        E_stored_initial= P_cap_storage/2;
    else
        E_stored_initial= sol.E_storage(end);
    end

    % Constraint: stored energy = eta*power_in - power_out/eta
    % NOTE: E_storage is stored energy at End of Hour
    storageBalance = optimconstr(num_hours);
    % for remaining hours, stored energy= previous stored energy + eta*power in - power out:
    storageBalance(1) = E_storage(1) == E_stored_initial + eta_storage*P_storage_in(1) - P_storage_out(1)/eta_storage;
    storageBalance(2:num_hours) = E_storage(2:num_hours) == E_storage(1:num_hours-1) +...
           eta_storage*P_storage_in(2:num_hours) - P_storage_out(2:num_hours)/eta_storage;
    linprob.Constraints.storageBalance = storageBalance;   

    % Constraint: power into storage over 1 hour cannot exceed free space in
    % storage at end of previous hour
 % Aug. 25, remove this constraint   
 %   cons_Pin_Ecap= optimconstr(num_hours);
 %   cons_Pin_Ecap(1:num_hours)= eta_storage*P_storage_in(1:num_hours) <= E_cap - E_storage(1:num_hours); 
 %   linprob.Constraints.PowerIn_Ecap= cons_Pin_Ecap;

    % Constraint: power out of storage cannot exceed stored energy
    cons_Pout_Estored= optimconstr(num_hours);
    cons_Pout_Estored(1)= P_storage_out(1)/eta_storage <= E_stored_initial;
    cons_Pout_Estored(2:num_hours)= P_storage_out(2:num_hours)/eta_storage <= E_storage(1:num_hours-1); 
    linprob.Constraints.PowerOut_Estored = cons_Pout_Estored;

    % Constraint: Power Balance, with curtailment
    cons_Pout_Demand= optimconstr(num_hours);
%     if P_cap_gas > 0
    cons_Pout_Demand(1:num_hours)= P_re - P_storage_in + P_storage_out + P_gas - P_curtail == P_demand;
%     else
%         cons_Pout_Demand(1:num_hours)= P_re - P_storage_in + P_storage_out - P_curtail == P_demand;
%     end
    linprob.Constraints.PowerOut_Demand = cons_Pout_Demand;

%      if P_cap_gas > 0
    % % Constraint: Gas plant has minimum power output when it is on
    cons_Gas_Pmin= optimconstr(num_hours);
    cons_Gas_Pmin(1:num_hours)= P_gas >= gas_on*P_min_gas;
    linprob.Constraints.gas_Pmin = cons_Gas_Pmin;

    % Constraint: Relate gas power ramping to gas power
    % gas power = gas_power + P_gas_ramp_up - power_out/eta
    % Note P_gas is power at end of hour
    gasRamp = optimconstr(num_hours-1);
    % for remaining hours, stored energy= previous stored energy + eta*power in - power out:
    %gasRamp(1) = E_storage(1) == E_stored_initial + eta_storage*P_storage_in(1) - P_storage_out(1)/eta_storage;
    gasRamp(2:num_hours) = P_gas(2:num_hours) == P_gas(1:num_hours-1) +...
        P_gas_ramp_up(2:num_hours) - P_gas_ramp_down(2:num_hours);
    linprob.Constraints.gasRamp = gasRamp;

    % https://math.stackexchange.com/questions/1851140/binary-integer-variables-in-linear-programming
    % Constraint: gas_on flag is 1 only when P_gas > 0
    % Algorithm tries to minimize cost, so it will set gas_on to 0 otherwise
    cons_Gas_On= optimconstr(num_hours);
    M= max(P_demand); % Just need large number
    cons_Gas_On(1:num_hours)= M*gas_on >= P_gas;
    linprob.Constraints.gas_On = cons_Gas_On;

    % z tracks the generator start-up
    w = optimexpr(num_hours); % Allocate w. Auxiliary variable
    % with spike at time gas turns on
    idx = 2:(num_hours);
    w(idx) = gas_on(idx) - gas_on(idx-1);
    w(1) = 0;
    switchcons = w - gas_turns_on <= 0;
    linprob.Constraints.switchcons = switchcons;
%      end
     
    %% Objective is to maximize profit
    % hourly profit = Hourly Revenue - Hourly cost
    % Hourly revenue= LMP * power provided to grid
    % power provided to grid= RE power - P sent to storage + P leaving storage

    % Cost_to_store is cost when power enters storage

%     Objec= sum(  LMP.*( P_re - P_storage_in + P_storage_out)...
%         -cost_to_store.*P_storage_in - cost_curtail.*P_curtail... % Van den Bergh uses LMP for Pcurtail
%         -cost_to_store.*E_storage...
%         -LMP.*P_gas - cost_ramp*(P_gas_ramp_up + P_gas_ramp_down)...
%         -cost_start*gas_turns_on );

    % JK Oct. 19, 2019. Just minimize entire electric supply costs
%     Objec= sum(  ...
%         -cost_to_store.*P_storage_in - cost_curtail.*P_curtail... % Van den Bergh uses LMP for Pcurtail
%         -cost_to_store.*E_storage...
%         -cost_ramp*(P_gas_ramp_up + P_gas_ramp_down)...
%         -cost_start*gas_turns_on );

    % JK Jan. 20, 2020
    Objec= sum(  ...
        -cost_to_store.*P_storage_in - cost_curtail.*P_curtail... % Van den Bergh uses LMP for Pcurtail
        ...% We decided we're not sure if there's a cost to keep energy stored -cost_to_store.*E_storage...
        -cost_ramp*(P_gas_ramp_up + P_gas_ramp_down)...
        -cost_start*gas_turns_on );
    
    

    linprob.Objective= Objec;
    %% Solution

    % Adjust options to speed up algorithm. Most of the default tolerances are
    % very small- a lot smaller than needed for our estimates
    % See: https://www.mathworks.com/help/optim/ug/intlinprog.html#btv2x05

    options = optimoptions('intlinprog',... 
        'ConstraintTolerance',  1e-3,...
        'CutMaxIterations',     40,...
        'MaxNodes',             500,...
        'MaxTime',              500,... % JK: set to essentially unlimited. error if returns empty soln when runs out of time
        'AbsoluteGapTolerance', 10,...
        'RelativeGapTolerance', .01,...
        'IntegerPreprocess',    'advanced',...
        'CutGeneration',        'advanced',...
        'RootLPAlgorithm',      'primal-simplex',...
        'HeuristicsMaxNodes',   100,...
        'Heuristics',           'advanced',...
        'ObjectiveImprovementThreshold',    1e-3,...
        'BranchRule',           'strongpscost'); %'fmincon');%;%, 'MaxIterations', 100*6*8761 );
    % options.Preprocess = 'none'; % workaround MATLAB Bug :-(

    [dispatchsol,annual_profit,exitflag,output] = solve(linprob, 'Options', options); %#ok<ASGLU>
    % Note: annual_profit here is the value of the objective function for control optimization.
    % It does not account for ICC and OPEX as used in LCOE calculation
    
    if exitflag ~= 1
        sol.exitflag= exitflag;
    end

    sol.P_storage_in=    [ sol.P_storage_in;     dispatchsol.P_storage_in ];
    sol.P_storage_out=   [ sol.P_storage_out;    dispatchsol.P_storage_out ];
    sol.E_storage=       [ sol.E_storage;        dispatchsol.E_storage ];
    sol.P_curtail=       [ sol.P_curtail;        dispatchsol.P_curtail ];
%     if P_cap_gas > 0
    sol.P_gas=           [ sol.P_gas;            dispatchsol.P_gas ];
    sol.gas_on=          [ sol.gas_on;           dispatchsol.gas_on ];         % Flag equals 1 when Gas supply is on
    sol.gas_turns_on=    [ sol.gas_turns_on;     dispatchsol.gas_turns_on ];   % Flag equals 1 when Gas turns on
    sol.P_gas_ramp_up=   [ sol.P_gas_ramp_up;    dispatchsol.P_gas_ramp_up ];  % For debugging
    sol.P_gas_ramp_down= [ sol.P_gas_ramp_down;  dispatchsol.P_gas_ramp_down]; % For debugging
%     else
%         % Dummy outputs so final results can show "0"'s for gas plant results
%         sol.P_gas= 0.*sol.P_storage_in;
%         sol.gas_on=  0.*sol.P_storage_in;
%         sol.gas_turns_on=  0.*sol.P_storage_in;
%         sol.P_gas_ramp_up=  0.*sol.P_storage_in;
%         sol.P_gas_ramp_down=  0.*sol.P_storage_in;
%     end
end

sol.E_cap= max(sol.E_storage);

% P_re= P_re_year';
% P_demand= P_demand_year;
% LMP= LMP_year;

% storage power losses at each hour              
sol.storage_losses=   (1-eta_storage)*sol.P_storage_in...      % P_loss enter storage
                + (1-eta_storage)/eta_storage* sol.P_storage_out;  % P_loss leave storage

% Power Supply-Demand balance is:
% P_demand == P_gas - P_storage_in + P_storage_out + P_re - P_curtail
% Required gas supply for power balance:
% P_gas= P_demand + P_storage_in - P_storage_out - P_re + P_curtail;

% NOTE: P_gas in the line above shows 

% When there is an oversupply of RE power and no storage space, 
% then some turbines/WECs needs to be shut down "curtailed"


if plotOn
    close all %#ok<UNRCH>
    figure;
    hours_to_plot= 1:730;
    hold on
    plot( P_re(hours_to_plot),'-+' );
    plot( P_demand(hours_to_plot),'-o' );
    plot( P_storage_in(hours_to_plot),'-s' );
    plot( P_storage_out(hours_to_plot),'-d' );
    plot( E_storage(hours_to_plot) );
    plot( E_cap + 0.*hours_to_plot );
    plot( P_curtail(hours_to_plot), 'r-x');
    plot(P_gas(hours_to_plot), '-o');
    plot(P_min_gas + 0.*hours_to_plot );
    plot(LMP(hours_to_plot), '-x');
    legend('RE Supply', 'Demand', 'Power into storage', 'Power out of storage',...
           'Stored Energy', 'Storage capacity', 'Curtailed Power', 'Gas Power', 'Min Gas Power', 'LMP');
end

end