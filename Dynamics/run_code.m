% Create response versus frequency plot for a single sea state.
% Create plot of response statistics
% For the unmodified and modified platform


%% 1. RAO's
% Select sea states for comparison
runType= 'RAO';
sea_state.Hvector= 1;           % significant wave height (m)
sea_state.wpVector= 2*pi/16;    % Peak frequency (rad/s)
sea_state.Uvector= 16;          % Mean wind speed (m/s)
sea_state.p= 1;                 % Fraction of life in sea state (-)
sea_state.seaType= 'Bret';      % Spectrum type
%%
% Baseline FWT
geoMod= {'Platform', '6MWPlatform'}; % Indicates modification wrt OC3 platform
run(geoMod, runType, sea_state);
%%
% Standalone RM3
geoMod= {'TMDe', 'Platform', 'RM3Platform', 'heavePlate'};
run(geoMod, runType, sea_state);
%%
% FWT with heave WEC
geoMod= {'TMDe', 'Platform', '6MWPlatform', 'heavePlate'};
run(geoMod, runType, sea_state);
%%
% Plot formatting
legend('Standalone FWT', 'Standalone WEC', 'Combined FWT-WEC');


%% lifetime response statistics
runType= 'iterate';
sea_state.Hvector= [1 1 1 3 3 3 6 6];               % significant wave height (m)
sea_state.wpVector= 2*pi./[8 11 16 8 11 16 11 16];  % Peak frequency (rad/s)
sea_state.Uvector= [8 8 8 16 16 16 20 20];          % Mean wind speed (m/s)
sea_state.p= [0.09 0.18 0.3 0.06 0.13 0.22 0.01 0.01];  % Fraction of life in sea state (-)
sea_state.seaType= 'Bret';                              % Spectrum type

% Baseline FWT
geoMod= {'Platform', '6MWPlatform'};
xss_baselineFWT= run(geoMod, runType, sea_state);

% Standalone RM3
geoMod= {'TMDe', 'Platform', 'RM3Platform', 'heavePlate'};
[xss_baselineWEC, FWT_WEC_relative_Velocity_baselineWEC]= run(geoMod, runType, sea_state);

% FWT with heave WEC
geoMod= {'TMDe', 'Platform', '6MWPlatform', 'heavePlate'};
[xss_combined, FWT_WEC_relative_Velocity_combined]= run(geoMod, runType, sea_state);

% Create plot
V3rel= [FWT_WEC_relative_Velocity_baselineWEC, FWT_WEC_relative_Velocity_combined];
X1=       [xss_baselineFWT(1) xss_combined(1)];
X3=       [xss_baselineFWT(3) xss_combined(3)];
X5=       [xss_baselineFWT(5) xss_combined(5)];
% Create plot:
annualStats(X1, X3, V3rel, X5)

pause= 1;
%%
function [xss, FWT_WEC_relative_Velocity] = run(geoMod, runType, sea_state)

%% INPUT MODIFICATIONS OF THE BASELINE PLATFORM
% geoMod options: 'TMDe', 'Platform', 'heavePlate', 'Sphere', 'WellsinSpar', 'none', 'OWC_Array', 'TMD' 'TMDeT'

%% WEC parameters - NREL RM3

rho= 1025;
g= 9.81;

externalFlag= 1;
mTMDe= 207.6e3; % WEC-Sim Verification paper
I55TMDe= 2.1e7; % WEC-Sim
%dJz= 2e6; %1.2e7*.25;%1.2e6; % Constant used in WEC-Sim example
dJVec= 1.2e6;% WEC-Sim
KJz= 0; % No heave stiffness coupling between WEC and FWT 

% WAMIT RM3 output. Assume wave period = infinity
A11TMDe= 2.292294E+02 * rho * externalFlag;
A15TMDe= 1.311723E+03 * rho * externalFlag;
A33TMDe= 1.984842E+03 * rho * externalFlag;
A51TMDe= 8.990518E+02 * rho * externalFlag;
A55TMDe= 2.369566E+04 * rho * externalFlag;

% Heave hydrostatic stiffness
outer_diameter= 20; % Neary paper
inner_diameter= 6;
area= pi * (outer_diameter^2 - inner_diameter^2) / 4;
C33TMDe= rho*g*area * externalFlag; % WEC-Sim rm3.out WAMIT data line 320

% Pitch hydrostatic stiffness
% This BM approximation assumes a deep structure
Itop_area= pi*( outer_diameter /2)^4/4;
Vsub= 2.5 * area;
BM= Itop_area/Vsub;
K= -2.5;
B= K/2; % center of buoyancy. Assume mid-keel

KB= B-K;
G= -.72; % RM3 verification paper
KG= G-K;
GM= BM - G + B;
C55TMDe= mTMDe*g*GM * externalFlag; %7347?

% Parameters for Surge Force: Input vectors of zS, xS, dVS, dAS
zSVector= -2.5/2; %mid-keel
xSVector= 0;
dVSVector= mTMDe/rho * externalFlag;
dASVector= A11TMDe * externalFlag;

% Heave Force: Input vectors of zH, xH, dVH, dAH <-- for strip theory OR Sphere
zHVector= -2.5; % keel
xHVector= 0;
dVHVector= mTMDe/rho * externalFlag;
dAHVector= A33TMDe * externalFlag;


% Coupling terms between FWT and WEC
dJx= 0;
dJy= 0;
% dJz= 1.2e6; % Constant used in WEC-Sim example
KJx= 0;
KJy= 0; % Constant used in WEC-Sim example
% KJz= 0;

% WEC parameters that are not used
kWells= 0;
VChamber= 0;
ICol= 0;
C33Col= 0;
ACol= 0;

rigidPitch= 1; %Equals 1 if WEC is rigidly connected to FWT in pitch

% WEC array placement
% angle= 0; %[0 1 2 3 4 5].*pi./3;
numTMDe= 1;
LxTMDe= 0; %-13.5.*cos(angle);
LyTMDe= 0; %13.5.*sin(angle);

dWellsx= 0;

%Location of junction point on FWT
Lfx= 0.*LxTMDe;
Lfy= 0.*LyTMDe;
Lfz= K/2;

%Location of WEC- may be in nacelle
Lwx= LxTMDe;
Lwy= LyTMDe;
Lwz= K/2;

xWellsForceSVector= zeros(1,numTMDe); %for water forcing hitting horizontal Wells
zWellsForceSVector= zeros(1,numTMDe); %for water forcing hitting Wells
xWellsForceHVector= zeros(1,numTMDe); %for water forcing hitting Wells
zWellsForceHVector= zeros(1,numTMDe); %For water forcing hitting bottom of column

controlKJxHeave= 0;
controlKJySurge= 0;
controlKJyHeave= 0;
controlKJzSurge= 0;


%% ITERATED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changing parameters:
%   For Surface plots: MVector is row. NVector is column.
%   For Line plots:  NVector is x-axis. MVector is in legend.

NVector= 1;%kWellsVector;
%NLabel= 'Wells turbine coefficient (Pa/m^3)';
MVector= dJVec; %AswlVector; %kznVector; %rFloatVector; % %%fVector; %dVector;
%MLabel= 'SWL area (m^2)'; %'Nonlinear spring coefficient, C_{33Non} (N/m^3)'; %%Nonlinear spring coefficient, C_{33Non} (N/m^3)'; %'Added SWL area (m^2)';
%   For Line plots:  NVector is x-axis. MVector is in legend.


%% COMPUTE OUTPUTS

for m= 1:length(MVector)   
    dJz= dJVec(m);
for n= 1:length(NVector)
 
% Baseline OC3 FWT platform parameters
%         Dfactor= 1;
%         Lfactor= 1;
%         Ltop= 14; %Keep the total length of the top segment at 14 m, where 4 m is submerged.
%         Lmid= Lfactor*8;
%         Lbot= Lfactor*108;
%         Dtop= 6.5;
%         Dbot= Dfactor*9.4;
%         tst= .0522; %For FWT. [m]

       if sum( strcmp( geoMod, 'heavePlate' ) )
            heavePlate= heavePlateParameters;
       else
           heavePlate= [];
       end
        
        %% Platform parameters
        if sum(strcmp(geoMod, '6MWPlatform'))
            % 5/25/2019
            % Define platform/tower parameters based on
            % parameters available in Beiter paper for 6 MW turbine and
            % NREL OC3 parameters for 5 MW spar.
            % Reduce platform depth to 90 m.

            % 6 MW turbine for Bieter paper
            D_rotor= 155; %[m]
            mEnd= 3.655e5; % [Kg]
            tower_height= 97; % [m]
            D_tower_top= 3.51;
            D_tower_bot= 6; %[m]
            D_SWL_platform= 7.8; %[m]
            Dtop= D_SWL_platform;
            D_bot_platform= 12; %[m]
            Dbot= D_bot_platform;
            thickness_platform= .0522; % [m]
            thickness_tower= .027; % [m]
            platform_depth= 90;

            % Try these values 5/28/2019
            % Code adjusts ballast to satisfy this submergence
            Ltop= 14;
            LtopSub= 10;
            Lmid= 8;
            Lbot= platform_depth - LtopSub - Lmid;

            % OC3
            %Ltop= 14;
            %Lmid= 8;
            %Lbot= 108;
            %LTank= 35.75; %Tank bottom is 4 m from bottom
            %Dtop= 6.5;
            %Dbot= 9.4;
            % tst= .0522; %[m]

            % Approximate Fthrust
            rho_a= 1.125;
            a= 1/3;
            Aswept= pi * (D_rotor/2)^2;
            U= 12; % wind speed (m/s)
            Fthrust= 0.5 * rho_a * Aswept * U^2 * 4*a*(1-a);

            mSphereSteel= 0;

        elseif sum(strcmp(geoMod, 'RM3Platform'))
            % RM3 platform parameters
            
            mEnd= 0; % [Kg] No turbine
            thickness_tower= 1e-20; % [m] No tower mass
            tower_height= 97; % [m]
            D_tower_top= 1e-10; % No tower - negligibly small
            D_tower_bot= 1e-10; % No tower - negligibly small
            
            D_SWL_platform= 6; %[m]
            Dtop= D_SWL_platform;
            D_bot_platform= 6; %[m]
            Dbot= D_bot_platform;
            thickness_platform= .0254; % [m] Neary PDF page 158
            platform_depth= 42;

            % Code adjusts ballast to satisfy this submergence
            Ltop= 2; % length above SWL [m]
            LtopSub= 5; % arbitrary length since constant diameter [m]
            Lmid= 8; % arbitrary length since constant diameter [m]
            Lbot= platform_depth - LtopSub - Lmid;
            
            mSphereSteel= 0;

        end
        
        %% Calculate platform parameters (like ballast) for neutral buoyancy
        
        % Calculate modification steel and concrete for cost. Neglect for experiments
        if sum( strcmp( geoMod, 'heavePlate' ) )
            mAdded_to_platform= heavePlate.mTMDe;
            z_of_mAdded_to_platform= heavePlate.Lwz; % platform depth
        else
            mAdded_to_platform= 0;
            z_of_mAdded_to_platform= 0;
        end
        
         [zCM, mPlatform, IplatformCMx, ITot_Swl_x, C55matplat, C33matplat,...
           mPlatformSteel, mPlatformConc, z_Platform, diameter_Platform,...
           zCMTower, mTowerSteel, z_Tower, diameter_Tower]...
            = platform_structure_and_hydrostatic_properties_modifyOC3spar...
            (D_bot_platform, D_SWL_platform, D_tower_bot, D_tower_top, mEnd,...
             Lbot, Lmid, LtopSub, Ltop,...
             tower_height, thickness_platform, thickness_tower, mAdded_to_platform, z_of_mAdded_to_platform);
         
         % Parameters for table in paper
         mPlatform_withPlate_tons= (mPlatform + mAdded_to_platform) / 1000 % {metric tons)
         mSteelPlatform_withPlate_tons= mPlatformSteel + mAdded_to_platform
         mConcretePlatform_tons = mPlatformConc / 1000
         z_cmPlatform_withPlate = (mPlatform_withPlate_tons*zCM + mAdded_to_platform*z_of_mAdded_to_platform)/...
             (mPlatform_withPlate_tons + mAdded_to_platform)
         
       %Pack
       platformParam= [zCM, LtopSub, mPlatform, IplatformCMx,...
           C55matplat, C33matplat, mPlatformSteel,...
           mPlatformConc, mSphereSteel, Ltop, Lmid, Lbot, Dtop, Dbot, 1]; %tst]; % 1 for OC
     
       %% WEC parameters
       
       TMDeParam= []; % initialize
       if any( strcmp( geoMod, 'TMDe') )
    %Geometry- wrt SWL at WEC's axis.
        TMDeParam.mTMDe= mTMDe;
    %     TMDeParam.mTMDe_surge= mTMDe_surge;
        TMDeParam.I55TMDe= I55TMDe; % Ixx
        % Iyy, Izz not excited in simulation.
        % Ixy= 0 symmetric about z axis.
        % Ixy= I45 will be 0 for symmetric array

        TMDeParam.A11TMDe= A11TMDe;
        TMDeParam.A15TMDe= A15TMDe;
        TMDeParam.A51TMDe= A51TMDe;
        TMDeParam.A55TMDe= A55TMDe;
        TMDeParam.A33TMDe= A33TMDe;

        TMDeParam.C33TMDe= C33TMDe;
        TMDeParam.C55TMDe= C55TMDe;

        % Location of junction point on FWT
        TMDeParam.Lfx= Lfx;
        TMDeParam.Lfy= Lfy;
        TMDeParam.Lfz= Lfz;

        % Location of WEC- may be in nacelle
        TMDeParam.Lwx= Lwx;
        TMDeParam.Lwy= Lwy;
        TMDeParam.Lwz= Lwz;

        %%


        %Wells parameters
        TMDeParam.dWellsx= dWellsx;

        %Pressure EOM
        TMDeParam.kWells= kWells;
        TMDeParam.VChamber= VChamber;

        %Water column EOM
        TMDeParam.ICol= ICol; %includes added mass
        TMDeParam.C33Col= C33Col;
        TMDeParam.ACol= ACol;


        %Parameters that may differ for each component of array
            %Each of these may be a horizontal vector
        TMDeParam.dJx= dJx;
        TMDeParam.dJy= dJy;
        TMDeParam.dJz= dJz;
        TMDeParam.KJx= KJx;
        TMDeParam.KJy= KJy;
        TMDeParam.KJz= KJz;
        etaJ= 0.85; % HARD-CODED INTO SIMULATEFWT_WEC.m % standard when include losses: .85;

        %For lateral and vertical forcing- depends on depth
        %2. For Surge Force: Input vectors of zS, xS, dVS, dAS
        TMDeParam.zSVector= zSVector;
        TMDeParam.xSVector= xSVector;
        TMDeParam.dVSVector= dVSVector;
        TMDeParam.dASVector= dASVector;
        %3. For Heave Force: Input vector of  zH, xH, dVH, dAH <-- for strip theory OR Sphere
        TMDeParam.zHVector= zHVector;
        TMDeParam.xHVector= xHVector;
        TMDeParam.dVHVector= dVHVector;
        TMDeParam.dAHVector= dAHVector;

        TMDeParam.xWellsForceSVector= xWellsForceSVector; %for water forcing hitting horizontal Wells
        TMDeParam.zWellsForceSVector= zWellsForceSVector; %for water forcing hitting Wells
        TMDeParam.xWellsForceHVector= xWellsForceHVector; %for water forcing hitting Wells
        TMDeParam.zWellsForceHVector= zWellsForceHVector; %For water forcing hitting bottom of column

        TMDeParam.controlKJxHeave= controlKJxHeave;
        TMDeParam.controlKJySurge= controlKJySurge;
        TMDeParam.controlKJyHeave= controlKJyHeave;
        TMDeParam.controlKJzSurge= controlKJzSurge;

        TMDeParam.kzn= 0; 
        TMDeParam.rigidPitch= rigidPitch; %Equals 1 if WEC is rigidly connected to FWT in pitch- for torque forcing
        TMDeParam.rigidTranslate= 1;% NOT USED %Equals 1 if want to override hinge matrix, using K1in

        % For Relative motion? 6/15/2019 not used
        TMDeParam.K1In= 0; % For rigid surge approximation K1In;  
        TMDeParam.K2In= 0; %1e15; %K2In;
        TMDeParam.K3In= 0; %1e15; %K3In;
        TMDeParam.d1In= 0; %d1In;
        TMDeParam.d2In= 0; %d2In;
        TMDeParam.d3In= 0; %d3In

        TMDeParam.dWellsz= 0; %dWellsz;
        TMDeParam.controlK11= 0; %controlK11; %overrides K1In
        TMDeParam.controlK33= 0; %controlK33; %Overrides K3In
        
        TMDeParam.etaJ= 0.8;
       end
       
   %% Calculate response
       
   iterateParams.geoMod= geoMod;
   iterateParams.cylParam= [];
   iterateParams.TMDeParam= TMDeParam;
   iterateParams.platformParam= platformParam;
   iterateParams.heavePlate= heavePlate;
   iterateParams.WellsParam= [];
   iterateParams.frameParam= [];
   iterateParams.nonLinParam= {0, 'lin', 'lin'};
   iterateParams.domainSwitch= 'frequency';
   iterateParams.runType= runType;
   iterateParams.sea_state= sea_state;
   
   
   % Powers are for entire array.
   [xss, xmax, DEL, DEL_stress, D_tot, maxStress, meanPower_kW, meanPower2_kW,...
       mean_V33_WEC_wrt_Platform_RMS,...
       maxPower_kW, windPowerOut, maxWindPower]=...
       iterate_windSeaStates_20190519(iterateParams); 

   % Statistic outputs
   
%     %DOFs: 10+3+TMD+2*TMDe+2*N
        X1out(m,n)= xss(1)
        X3out(m,n)= xss(3)
        X5out(m,n)= xss(5) * 180/pi % [deg]
        FWT_WEC_relative_Velocity(m,n)= mean_V33_WEC_wrt_Platform_RMS; % [m/s]
        Powerout(m,n)= meanPower_kW %[kW]

end

end

end


function heavePlate= heavePlateParameters

    rho= 1025;

    heavePlate.mTMDe= 244.7e3; % (kg)

    R= 15; % Disk radius [m]

    % Standard formula for thin disk mass moment of inertia:
    heavePlate.I55TMDe= 0.25 * heavePlate.mTMDe * R^2; % (kg)

    heavePlate.Lwz = -90; % Location wrt SWL (m)
    heavePlate.Lfz= -90; % Location of junction point on FWT
    
    % Heave plate added mass parameters are based on Newman text page 145.
    % Added mass parameters are wrt plate centroid:
    heavePlate.A11TMDe= 0;
    heavePlate.A33TMDe= pi * rho * R^2;
    heavePlate.A55TMDe= 0.125 * pi * rho * R^2;

end







