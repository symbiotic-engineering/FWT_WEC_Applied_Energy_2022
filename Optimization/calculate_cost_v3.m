function [ ICC, LCOE ]= calculate_cost_v3(n_turbines, n_WECs, storage_E_cap, storage_Power_cap, mean_RE, geology, AEP_one_WEC, AEP_one_WT)
% calculate capital cost and LCOE for wind-wave-storage farm, including
% offshore transmission system.
% Inputs are:
% n_turbines: number of 6 MW floating turbines
% n_WEC's: number of 286 kW WEC's.
    % Assume that n_turbine of the n_WEC's are attached to turbines
% storage_E_cap: capacity of energy storage system [MWh]
% max_RE: maximum power output by RE system BEFORE transmission losses [MW]
    % (used for calculating transmission cost)
% mean_RE: mean power output by RE system AFTER all losses
    % (used for calculating LCOE)
% geology: either 'aquifer' or 'rock'. 'aquifer' is much cheaper.
 
% Update by MNH on 9/3/2020 to change costs to 2020 USD

distance_to_shore = 50; % km - MNH added 11/2/2020 for electrical system cost calculation
euro_2014_to_usd_2014 = 1/0.784;

% for storage costs
inflation_2009_to_2020 = 1.1426;
inflation_2009_to_2020 = 1.2077;

% for WEC costs
inflation_2014_to_2020 = 1.0354; 
inflation_2014_to_2018 = 1.0607;
inflation_2014_to_2020 = 1.0945;

% for WT costs
inflation_2015_to_2020 = 1.0342; 
inflation_2015_to_2018 = 1.0594;
inflation_2015_to_2020 = 1.0932;

FCR = 11.4/100; % fixed charge rate %
distance_to_shore = 50; % [km]

n_standalone_WECs= max(0, n_WECs - n_turbines);
n_attached_WECs= min(n_turbines, n_WECs);
n_standalone_turbines = max(0, n_turbines - n_WECs);

P_cap_WEC= 286; % [kW]
P_cap_WT = 6; % [MW]

WEC_cf = 0.30; % WEC capacity factor = 30%
eta1_WEC = 0.95; % WEC efficiency based on losses due to device availability = 95%
% eta2_WEC = 0.98; % WEC efficiency based on transmission losses = 98% -
% NOT USED because transmission losses computed
eta2_WEC = 1;

WT_cf = 0.30; % wind turbine capacity factor = 30%
WT_uptime = 0.95; % Wind turbine uptime = 95%

%% Stand alone WEC costs - MNH Nov 23, 2020 now based on 100 unit costs not exponential curves
if n_standalone_WECs ~= 0
    % NH Nov 23, 2020 now based on 100 unit costs not exponential curves
    ICC_WEC_cost_per_kw = 13630*inflation_2014_to_2020; % [2020 $] 
    AOE_WEC_cost_per_kw = 330*inflation_2014_to_2020; % [2020 $/yr]
    
    % OLD when was based on exponential cost curves
    % P_cap_WEC must be in kW for this equation, P_cap_WEC*n_standalone_WECs = nameplate capacity (kW)
    % ICC_WEC_cost_per_kw = 396654*(P_cap_WEC*n_standalone_WECs)^(-0.332); % [2020 $] 
    % AOE_WEC_cost_per_kw = 99585*(P_cap_WEC*n_standalone_WECs)^(-0.55); % [2020 $/yr] 
    
    % JK July 27, 2019: 
    ICC_WEC= ICC_WEC_cost_per_kw * n_WECs * P_cap_WEC;
    AOE_WEC= AOE_WEC_cost_per_kw * n_WECs * P_cap_WEC;
    
else
    ICC_WEC = 0;
    AOE_WEC = 0;
end

AEP_WEC = P_cap_WEC*n_standalone_WECs*WEC_cf*eta1_WEC*eta2_WEC*365*24/10^3; % [MWh]

%% Stand alone WT costs
if n_standalone_turbines ~=0
    % MH 11/2/2020 - this is not matching the excel sheet. Updated below. ICC_WT = 35272665*inflation_2015_to_2020*n_standalone_turbines; % [2020 $] from Wind LCOE excel sheet
    LostRev_WT = 150*AEP_one_WT;
    ElecCapEx_WT = (4.6022*distance_to_shore+208.61)*10^6/100-LostRev_WT;
    ICC_WT =  (29275666+ElecCapEx_WT)*inflation_2015_to_2020*n_standalone_turbines; % [2020 $] from Wind LCOE excel sheet
    AOE_WT = (4.2167*log(distance_to_shore)+122.54)/100*10^6*inflation_2015_to_2020*n_standalone_turbines; % [2020 $/yr]
else
    ICC_WT = 0;
    AOE_WT = 0;
end

AEP_WT = P_cap_WT*n_standalone_turbines*365*24*WT_cf*WT_uptime; % [MWh]

%% Combined WEC and WT costs - MNH Nov 23, 2020 now based on 100 unit costs not exponential curves
if n_attached_WECs ~=0 
    % NH Nov 23, 2020 now based on 100 unit costs not exponential curves
    ICC_WW_WEC_cost_per_kw = 7550*inflation_2014_to_2020; % [2020 $] 
    AOE_WW_WEC_cost_per_kw = 300*inflation_2014_to_2020; % [2020 $/yr]
    
    % OLD when was based on exponential cost curves
    % P_cap_WEC must be in kW for this equation
    % ICC_WW_WEC_cost_per_kw = 154699*(P_cap_WEC*n_attached_WECs)^(-0.295); % [2020 $] already in 2020 from Excel calcs
    % AOE_WW_WEC_cost_per_kw = 79015*(P_cap_WEC*n_attached_WECs)^(-0.538); %[2020 $/yr] already in 2020 from Excel calcs
    
    % JK July 27, 2019:
    %ICC_WW_WEC= ICC_WW_WEC_cost_per_kw * n_WECs * P_cap_WEC;
    %AOE_WW_WEC= AOE_WW_WEC_cost_per_kw * n_WECs * P_cap_WEC;
    % MNH Nov 2, 2020: should be "n_attached_WECs" not all of n_WECs
    ICC_WW_WEC= ICC_WW_WEC_cost_per_kw * n_attached_WECs * P_cap_WEC;
    AOE_WW_WEC= AOE_WW_WEC_cost_per_kw * n_attached_WECs * P_cap_WEC;
    
    % MH 11/2/2020 - this is not matching the excel sheet. Updated below. ICC_WW_WT = 35351601*inflation_2015_to_2020*n_attached_WECs; % [2020 $]
    LostRev_WW = 150*(AEP_one_WT+AEP_one_WEC);
    ElecCapEx_WW = (4.6022*distance_to_shore+208.61)*10^6/100-LostRev_WW;
    ICC_WW_WT =  (29343702+ElecCapEx_WW)*inflation_2015_to_2020*n_attached_WECs; % [2020 $] From Wind+WEC_LCOE excel sheet
    AOE_WW_WT = (4.2167*log(distance_to_shore)+122.54)/100*10^6*inflation_2015_to_2020*n_attached_WECs; %[2020 $/yr]
    
    AEP_WW_WEC = P_cap_WEC*n_attached_WECs*WEC_cf*eta1_WEC*eta2_WEC*365*24/10^3; % [MWh]
    AEP_WW_WT = P_cap_WT*n_attached_WECs*365*24*WT_cf*WT_uptime; % [MWh]
    
    ICC_WW = ICC_WW_WEC + ICC_WW_WT;
    AOE_WW = AOE_WW_WEC + AOE_WW_WT;
    
    AEP_WW = AEP_WW_WEC + AEP_WW_WT;
else
    ICC_WW = 0;
    AOE_WW = 0;
end
%% Installed Capital Cost


% Storage capital cost
if strcmp(geology, 'rock')
    % Convert  from [MWh] to [kWh]
    % Zakeri, total capital cost = 893 euros/kW average
    ICC_storage= storage_E_cap * 1000 * 30 * inflation_2009_to_2020 ... % $30/kWh for rock
        + 893 * euro_2014_to_usd_2014 * inflation_2014_to_2020 * 1000 * storage_Power_cap; % $1.2M/MW
else % aquifer
    ICC_storage= storage_E_cap * 1000 * 0.1 * inflation_2009_to_2020 ... % $0.1/kWh for aquifer
        + 893 * euro_2014_to_usd_2014 * inflation_2014_to_2020 * 1000 * storage_Power_cap; % $1.2M/MW
end

% Storage operating cost
AOE_storage = 2*inflation_2015_to_2020*storage_E_cap; % $2/MWh from Liu, in 2015 $, supported by Zakeri

% Total installed capital cost
ICC= ICC_WEC + ICC_WT + ICC_WW + ICC_storage;

% Total annual operating cost
AOE = AOE_WEC + AOE_WT + AOE_WW + AOE_storage;

% annual energy production
AEP = mean_RE*8766; % [MWh]

LCOE= ( ICC * FCR + AOE ) / AEP; %[ $ / MWh ]

end



