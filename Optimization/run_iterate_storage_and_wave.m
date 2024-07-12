function sol= run_iterate_storage_and_wave(varargin)
% RUN_ITERATE_STORAGE_WAVE calculates the wind-wave-storage farm / gas plant
% optimized time series and performance statistics of power supply, demand,
% losses, and cost.
%
% The statistics are normalized in terms of peak power or energy.
% Peak demand is scaled to equal peak wind/wave supply.
%
% The energy storage capacity resulting from optimization is an output.
%
% Inputs: N_WEC_VEC, P_CAPACITY_STORAGE, RE_PENETRATION.
% N_WEC_VEC, P_CAPACITY_STORAGE may be scalars or vectors.
% The # of 6 MW wind turbines is fixed to 100.
%
% If RUN_ITERATE_STORAGE_WAVE is called without arguments, then the inputs
% listed in lines 28-33 are used.
% plots are only generated when RUN_ITERATE_STORAGE_WAVE is called without
% arguments.


tic; % time the run.

%%%%% ------ Inputs ------ %%%%%
if nargin > 0
   % Fraction of renewable energy from Wave. Current fixed parameters try
   % to adjust n_WEC and n_turbines so P_farm is 600 MW.
   frac_RE_wec= varargin{1};
   P_Capacity_Storage_frac_vec= varargin{2};
   RE_penetration= varargin{3};
   generatePlots= 0;
else
    % n_WEC_vec, number of 0.286 MW WEC's. (n_turbines of the n_WEC's are attached to turbines. the remainder are standalone).
    % n_turbines, number of 6 MW turbines  
    
    % Fraction of renewable energy from Wave. Current fixed parameters try
    % to adjust n_WEC and n_turbines so P_farm is 600 MW.
    frac_RE_wec= [0.5 1]; %[0:0.1:1];
    % n_WEC_vec= [0:300:3000]; OLD 

    P_Capacity_Storage_frac_vec= [0 1]; %[0 0.05 0.1 0.2:0.2:1]; % fraction of Community peak demand 
    % E_Capacity_hrs_vec= [0 1 2 3 4 5 10 24 36 48]; %linspace(10, 100, 10); % Storage Capacity in [# FWT-WEC hourly energy capacity] 
    RE_penetration= 0.5; % fraction of grid supplied by renewable energy (wind+ wave)
    generatePlots= 1;
end

%% Fixed parameters in our analysis
    
% number of 6 MW wind turbines for nameplate capacity of 600 MW farm
% n_turbines= 100; 
P_cap_farm_target= 600; % [MW]

dist_to_shore = 50; % distance of farm to shore [km]
transmission = 'HVAC'; % If HVDC or HVAC

P_cap_turbine= 6; % [MW]
P_cap_WEC= 0.286; % [MW]
array_soil_losses= .085; % Array and Soiling losses
%   JK 8/16/2019: we have been using 12.5%, but Fingersh paper says just 8.5% losses

geology= 'aquifer'; % geology for storage system. either 'aquifer' or 'rock'. 'aquifer' is much cheaper.

eta_storage= 0.9; % efficiency when power enters and exits storage
cost_to_store= 2; % For now, just assume $2/MWh when entering storage, sources cited in manuscript


%% %%% ------ Derived parameters ------ %%%%%

% vector of fraction of RE capacity that is from waves (the rest is from wind). Used in final output plot
% frac_RE_wec= n_WEC_vec .* P_cap_WEC ./ (n_WEC_vec .* P_cap_WEC + n_turbines * P_cap_turbine );


%% Iterate through n_WEC_vec and E_capacity_hrs, and calculate results

% Initialize Matrices of final outputs
%Average system efficiency:
    %Annual Energy produced by storage system /Annual wind power
E_efficiency= [];
%maximum adjustment that plant must make in 1 hour relative to hour before
max_frac_gas_diff_per_hour= []; 

%minimum adjustment that plant must make in 1 hour relative to hour before
min_frac_gas_diff_per_hour= []; 
    
for i= 1:length(P_Capacity_Storage_frac_vec) %(E_Capacity_hrs_vec) 
    P_Cap_Store_Frac= P_Capacity_Storage_frac_vec(i);
    
    for j= 1:length(frac_RE_wec) 
        frac_RE_wec_j= frac_RE_wec(j);
        clc
        i
        j %#ok<*NOPRT>
        %n_WECs= n_WEC_vec(j);
        
        % Solve for # WEC's and # turbines for P_farm = 
        
        % Ax= b
        A= [P_cap_WEC                     P_cap_turbine
            P_cap_WEC*(frac_RE_wec_j -1)  P_cap_turbine*frac_RE_wec_j];
        b= [P_cap_farm_target; 0];
        x= A\b;
        n_WECs= round(x(1)); % Always at least 1 WEC
        n_turbines= round(x(2)); % Always at least 1 Wind turbine
        
        P_cap_wind_farm= n_turbines * P_cap_turbine; % [MW]
        P_cap_wave_farm= n_WECs * P_cap_WEC; % [MW]
        P_cap_farm= P_cap_wind_farm + P_cap_wave_farm;
        
        % E_cap now a variable. E_cap= E_Capacity_hrs_vec(i)*( P_cap_farm ); % [MWh]
        % Assumption: Limit P_storage_in and P_storage_out to be P_farm
        
        % Optimize the control of the system.
        sol = get_dynamics_performance(n_turbines, n_WECs, P_cap_farm, P_Cap_Store_Frac, RE_penetration,...
        cost_to_store, eta_storage, dist_to_shore, transmission, array_soil_losses);
    
        sol.mean_Demand= mean(sol.Demand);
        sol.RE_energy= sum(sol.RE_Supply_transmitted);
        
        P_mean_storage_losses(j,i)= mean( sol.P_storage_in - sol.P_storage_out );
        P_mean_curtail(j,i)= mean( sol.P_curtail );
        % P_final_losses is all losses besides transmission losses.
        % mean_Re already accounts for transmission losses.
        P_final_losses= sol.P_storage_in - sol.P_storage_out + sol.P_curtail;
        P_mean_RE(j,i)= sol.mean_RE; 
        
        % sol.mean_RE is mean RE transmitted AFTER transmission losses
        % June 2021. Farm_CF(j,i)= sol.mean_RE./P_cap_farm; %#ok<*AGROW>
        Farm_CF(j,i)= (sol.mean_RE - mean(P_final_losses))./P_cap_farm;
        Storage_Out_CF(j,i)= sol.storage_power_mean./sol.storage_power_cap;
        Storage_E_capacity(j,i)= sol.E_cap./sol.Demand_peak; %./P_cap_farm; 
        Storage_P_capacity(j,i)= sol.storage_power_cap./sol.mean_Demand;
        Storage_P_mean(j,i)= sol.storage_power_mean./sol.mean_Demand;
        E_efficiency(j,i)= sol.E_efficiency;
        RE_capacty(j,i)= sol.max_RE; % RE power capacity BEFORE losses
        Gas_capacity(j,i)= sol.gas_capacity./sol.mean_Demand; % Gas capacity is increased as needed to meet demand
        max_frac_gas_diff_per_hour(j,i)= sol.max_frac_gas_diff_per_hour;
        min_frac_gas_diff_per_hour(j,i)= sol.min_frac_gas_diff_per_hour;
        avg_frac_gas_diff_per_hour(j,i)= sol.avg_gas_ramp;
        std_frac_gas_diff_per_hour(j,i)= sol.std_gas;
        gas_power_fluctuation_ratio(j,i)= (max(sol.P_gas) - min(sol.P_gas))/mean(sol.P_gas);
        min_gas_frac_of_cap(j,i)= sol.min_gas_frac;
        gas_CF(j,i)= sol.gas_CF;
        E_curtail(j,i)= sol.E_curtail./sol.RE_energy;
        exitFlag(j,i)= sol.exitflag; % track how solver exits.
        
        AEP_one_WEC= sum(sol.Wave_Supply_transmitted)/n_WECs; % [MWh]
        AEP_one_WT= sum(sol.Wind_Supply_transmitted)/n_turbines; % [MWh]
        
        if isnan(AEP_one_WEC) || isinf(AEP_one_WEC)
            AEP_one_WEC= 0;
        end
        if isnan(AEP_one_WT) || isinf(AEP_one_WT)
            AEP_one_WT= 0;
        end
        
        P_mean_curtailed= mean( sol.P_curtail );
        
        % Calculate costs - MNH updated Nov 23, 2020 to use v3 of cost calculation code, no longer exponential WEC costs, instead based on 100 unit costs
            [ capital_cost, lcoe ]= calculate_cost_v3(n_turbines, n_WECs, sol.E_cap, sol.storage_power_cap, sol.mean_RE-P_mean_curtailed, geology, AEP_one_WEC, AEP_one_WT);
        LCOE_RE(j,i)= lcoe; % [$/MWh]
        Cap_cost_RE(j,i)= capital_cost./P_cap_farm; % [$/MW]
    end

end

%% Plot Overall performance
close all
if  generatePlots
    figure;
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, Farm_CF, {'\fontsize{9}Wind-Wave-Storage Farm'; 'Capacity Factor after Losses'}, 3);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, Cap_cost_RE, {'\fontsize{9}Wind-Wave-Storage Farm'; 'ICC ($/MW)'}, 2);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, LCOE_RE, {'\fontsize{9}Wind-Wave-Storage Farm'; 'LCOE ($/MWh)'}, 1);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, Storage_Out_CF, {'\fontsize{9} Storage Output Capacity Factor'}, 6);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, Storage_E_capacity, {'\fontsize{9} Storage Energy Capacity'; '(average demand hours)'}, 5);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, E_efficiency, {'\fontsize{9}Wind-Wave-Storage Farm'; 'Electricity-to-grid efficiency'}, 4);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, gas_CF, '\fontsize{9}Gas Plant Capacity Factor', 9);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, gas_power_fluctuation_ratio, {'\fontsize{9} Gas Power Supply Fluctuation Ratio'}, 8);
    plot_final(P_Capacity_Storage_frac_vec, frac_RE_wec, E_curtail, {'\fontsize{9} Annual Curtailed Energy'; '/ Renewable Energy Produced'}, 7);
end

run_time= toc

%% save data
time= datestr(now, 'yyyy_mm_dd_HH_MM');
filename= ['Preliminary_Results\results_' time];
save(filename)

end

%Plot System Cost versus Energy storage capacity:
%wind farm- Navigant
%transmission line- Bertsekas-Negra
%converter stations- Bersekas-Negra
%Storage system- assume compressed air- aquaduct- Carnegie

%% ----------- HELPER FUNCTIONS ----------
function plot_final(E_Capacity_hrs, frac_RE_wec, data, title_str, n)
    subplot(3, 3, n);
    [x,y]= meshgrid(E_Capacity_hrs, frac_RE_wec);
    surf(x, y, data, 'FaceColor', 'interp', 'LineStyle', 'none');
    view([0 90])

    if n >= 7
        %xlabel('\fontsize{9}Storage capacity hours (hr)');
        xlabel({'\fontsize{9}Storage power capacity', '\fontsize{9}(fraction of peak demand)'});
    end
    if n==1 || n==4|| n==7
        ylabel({'\fontsize{9}Fraction of Renewable'; 'Energy from Wave'});
    end
    title(title_str);
    axis([min(E_Capacity_hrs) max(E_Capacity_hrs) min(frac_RE_wec) max(frac_RE_wec)]);
    colorbar
end
    
