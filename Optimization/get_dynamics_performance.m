function sol= get_dynamics_performance...
    (n_turbines, n_WECs, P_farm, P_Cap_Store_Frac, fraction_RE,...
    cost_to_store, eta_storage, dist_to_shore, transmission, soil_losses)       
% OUTPUTS: -- all stored in the struct "sol" --
%   Performance Statistics:
%       max_RE, mean_RE
%       gas_capacity, max_frac_gas_diff_per_hour, min_frac_gas_diff_per_hour,
%       min_gas_frac, avg_gas_ramp, std_gas, gas_CF,...
%       E_efficiency, E_curtail,...
%       E_cap, storage_power_cap, storage_power_mean,... 
%       Demand_peak, Demand_annual_energy.
%   Optimized time series:
%       P_curtail, P_gas, P_storage_in, P_storage_out, E_storage, storage_losses, E_cap

%%
plotOn= 0; %Set to 1 to make intermediate plots, set to 0 to not plot

% Lines below may be commented out and resulting vectors loaded below to reduce runtime
%%%%%%%%%%%%%%%%%%%%%%% GET FARM SUPPLY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Supply and Demand data scaled so that peak supply and peak demand both =
         %500.8 MW (Power output of 159 5 MW wind turbines with 21% losses)
         
         
[date_S, Wind_Supply, windSpeed, Wave_Supply]= power_supply_fun(n_turbines, n_WECs, soil_losses);

if plotOn
    plot_correlations_month_resolution(date_S, Wind_Supply, Wave_Supply); 
    plot_correlations_day_resolution(date_S, Wind_Supply, Wave_Supply);
end

% interpolate data to match every hour since data seems to be at slightly different times
date= datenum('1-Jan-2017'):1/24:datenum('31-Dec-2017'); %each integer is 1 day
% date = 1 data point per hour so all data is hourly
Wind_Supply= interp1(date_S, Wind_Supply, date, 'linear', 'extrap');
Wave_Supply= interp1(date_S, Wave_Supply, date, 'linear', 'extrap');

RE_Supply= Wind_Supply + Wave_Supply;
max_RE= max(RE_Supply); %RE capacity BEFORE losses
wind= interp1(date_S, windSpeed, date, 'linear', 'extrap'); %#ok<NASGU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIBE LOSSES - BEFORE STORAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(transmission, 'HVDC')
    % FOR HVDC
    % transmission line losses
    % Table of [Wind power before offshore substation, power delivered to grid...
    % P_into CS, power lost by CS, Power into cable, power lost by cable]
    [~, ~, P_in_CS, loss_CS, P_in_cab, P_loss_cab]= VSC(max_RE);
    
    %Wind farm provides Supply to CS (already accounts for 12.5% turine losses)
    loss_CS1= interp1(P_in_CS, loss_CS, RE_Supply); %power lost by CS
    P_out_CS1= RE_Supply-loss_CS1;
    P_out_CS1(P_out_CS1<0)= 0; %numerical error
    loss_cable= interp1(P_in_cab, P_loss_cab, P_out_CS1);
    P_out_cable= P_out_CS1-loss_cable;
    P_out_cable(P_out_cable<0)= 0; %numerical error
    loss_CS2=  interp1(P_in_CS, loss_CS, P_out_cable);
    P_out_CS2= P_out_cable - loss_CS2;
    P_out_CS2(P_out_CS2<0)= 0; %#ok<NASGU> %numerical error
    
    all_transmission_losses= loss_CS1 + loss_cable + loss_CS2;
else
    % FOR HVAC:
    % Total transmission losses from Todorovic pg 122 for 500 MW farm, 3x
    % 220kV, 50km from shore using HVAC (includes cable and transformer losses)
    P_farm_Capacity_round = round(P_farm,-2); % round P_farm to the nearest 100 MW so it works with code below
    P_loss_perc= HVAC_trans_losses(P_farm_Capacity_round, RE_Supply, windSpeed, dist_to_shore);
    
    all_transmission_losses = P_loss_perc.*RE_Supply;
    Wind_Supply_transmitted= (1-P_loss_perc).*Wind_Supply; % For plotting
    Wave_Supply_transmitted= (1-P_loss_perc).*Wave_Supply; % For plotting
end
RE_Supply_transmitted= (RE_Supply - all_transmission_losses)';
mean_RE= mean(RE_Supply_transmitted); % mean renewable energy AFTER losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%% GET POWER DEMAND DATA %%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Demand power is scaled by:
% E_Supply_RE = fraction_RE * E_Demand
% annual energy produced does not account for storage losses  

% unscaled demand power
unscaled_demand= demand_fun(date_S);
unscaled_demand= interp1(date_S, unscaled_demand, date,'linear', 'extrap');

E_Supply_RE= sum(RE_Supply_transmitted); % MWhr
E_Demand_required= E_Supply_RE /fraction_RE; % MWhr
E_Demand_unscaled= sum(unscaled_demand); % MWhr
Demand_scaling_factor= E_Demand_required / E_Demand_unscaled;

% Scaled demand power
Demand = Demand_scaling_factor * unscaled_demand;

Demand_peak= max(Demand);
Demand_annual_energy= sum(Demand);

%% SCALE STORAGE POWER CAPACITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_cap_storage= P_Cap_Store_Frac * Demand_peak;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotOn
    plot_correlations_month_and_day_resolution(date_S, Wind_Supply, Wave_Supply, Demand); %#ok<*UNRCH>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotOn==1
    %plot wind speed and power output 
    tickDates= datenum( 2017.*ones(3,1), [1 5 9]', ones(3,1) );
    DateString= datestr(tickDates, 'mmm ''yy');
    figure;
    plot(date, Wind_Supply,date, Wave_Supply, date, Demand); %power from farm into CS offshore
    title('Power Transmitted and Demand');
    ylabel('Power (MW)');
    set(gca, 'XTick', tickDates);
    set(gca, 'XTickLabel', DateString);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Load data instead of computing
% load historical_data5.mat date Demand Supply


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAS PLANT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%base load of gas plant

num_hrs_in_yr= 24*365;
base_gas= (1-fraction_RE) * E_Demand_required/ num_hrs_in_yr; % expected average power

CF_gas= 0.5; 
% standard value https://www.eia.gov/todayinenergy/detail.php?id=25652

% if fraction_RE == 1
%     P_cap_gas= 0;
% else
P_cap_gas= max(Demand); %base_gas/CF_gas;
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT TRANSMITTED SUPPLY AND DEMAND- NO STORAGE %%%%%%%%%%%%%%%%%%%%%%%%%
tickDates= datenum([[2010.*ones(3,1), [1 5 9]', ones(3,1)]; [2011.*ones(3,1), [1 5 9]', ones(3,1)]; [2012.*ones(3,1), [1 5 9]', ones(3,1)]]);
DateString= datestr(tickDates, 'mmm ''yy');

if plotOn== 1

    %Plot transmitted supply and demand- No Storage to show mismatch
    figure(11);
    plot(date, base_gas.*ones(1,length(date)), 'g', date, RE_Supply_transmitted+base_gas, 'b', date, Demand, 'r')
                %supply_transmitted = supplied wind power minus all losses NOT considering any storage losses
    legend('Gas Base Load', 'Transmitted Wind Power + Gas Base Load', 'Consumer Demand');
    hold on
    H=area(date, base_gas.*ones(1,length(date)));
    h=get(H,'children');
    set(h,'FaceColor','green', 'LineStyle', 'none'); %#Tada!
    plot(date, Demand, 'r');
    plot(date, base_gas.*ones(1,length(date)), 'g', 'LineWidth', 5);
    plot(date, RE_Supply_transmitted+base_gas, 'b');
    plot(date, Demand, 'r');

    xlabel('Date');
    ylabel('Power (MW)');
    title('Wind farm power and consumer demand- No energy storage');
    set(gca, 'XTick', tickDates);
    set(gca, 'XTickLabel', DateString);
    axis([date(1) date(end) 0 1.1*max(Demand)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate energy storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[date_LMP, LMP_data]= get_hourly_Electricity_Prices;

LMP= interp1(date_LMP, LMP_data, date, 'linear', 'extrap');

% Determine P_curtail, P_gas, P_storage_in, P_storage_out, E_storage, storage_losses, E_cap.
% This data is stored in the outputted data struct, sol.
sol= optimize_storage_control_for_profit(Demand, RE_Supply_transmitted, P_cap_gas, P_cap_storage, cost_to_store, LMP, eta_storage);

% Other outputs
sol.max_RE= max_RE;
sol.mean_RE= mean_RE;
sol.Demand_peak= Demand_peak;
sol.Demand_annual_energy= Demand_annual_energy;
sol.Wind_Supply_transmitted= Wind_Supply_transmitted;
sol.Wave_Supply_transmitted= Wave_Supply_transmitted;
sol.RE_Supply_transmitted= Wind_Supply_transmitted + Wave_Supply_transmitted;
sol.Demand= Demand;

% loss_Store is Power losses of storage system

total_losses= all_transmission_losses' + sol.storage_losses;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS AND OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Scalar outputs
    E_farm_supply= sum(RE_Supply(~isnan(RE_Supply))); % total RE supply before losses
    E_lost= sum(total_losses( ~isnan(RE_Supply) ) ); % all lost energy over entire year
    sol.E_efficiency= (E_farm_supply-E_lost)/E_farm_supply;

    % CEEP= mean(abs(P_gas - base_gas)); %mean adjustment required to gas plant
    diffCEEP= diff(sol.P_gas - base_gas, 1); %production difference between consecutive hours

    sol.gas_capacity= max(sol.P_gas);
    mean_gas= mean(sol.P_gas);
    sol.max_frac_gas_diff_per_hour= max(abs(diffCEEP))/sol.gas_capacity;
    sol.min_frac_gas_diff_per_hour= min(abs(diffCEEP))/sol.gas_capacity;
    sol.gas_CF= mean_gas/sol.gas_capacity;
    sol.avg_gas_ramp= mean( abs( diff(sol.P_gas) ) )/sol.gas_capacity;
    sol.std_gas= std( sol.P_gas ) / sol.gas_capacity;
    sol.storage_power_cap= max( max(sol.P_storage_in), max(sol.P_storage_out) ); % should be P_farm
    sol.storage_power_mean= mean( sol.P_storage_out ); % Just look at power out of storage
        % to indicate how much power flows through the storage system
    
    sol.min_gas_frac= min(sol.P_gas)/sol.gas_capacity;

    sol.E_curtail= sum( sol.P_curtail );

end





