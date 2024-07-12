function [Xssmean, Xssmax, DEL, DEL_stress, D_tot, maxStress, meanPower, meanPower2, mean_V33_WEC_wrt_Platform_RMS, ...
    maxPower_KW, windPowerOut, maxWindPower]= iterate_windSeaStates_20190519(varargin)
% Compute hinge parameters

%% POSSIBLE INPUTS
%if this function is called from command window, ensure that geoMod is 'none'
if isstruct(varargin{1})
    iterateParams= varargin{1};
    geoMod= iterateParams.geoMod;
    cylParam= iterateParams.cylParam;
    TMDeParam= iterateParams.TMDeParam;
    platformParam= iterateParams.platformParam;
    heavePlate= iterateParams.heavePlate;
    frameParam= iterateParams.frameParam;
    
    nonLinParam= iterateParams.nonLinParam;
        Cd= nonLinParam{1};
        mSwitch= nonLinParam{2};
        cSwitch= nonLinParam{3};
    domainSwitch= iterateParams.domainSwitch;
    runType= iterateParams.runType;
    
    useWpAB= 0;
    
    if ~isempty(iterateParams.WellsParam)
        WellsParam= iterateParams.WellsParam;
    else
        WellsParam= [0 0];
    end
    
    if ~isempty(iterateParams.sea_state)
        SeaStates= iterateParams.sea_state;
        Hvector= SeaStates.Hvector;
        wpVector= SeaStates.wpVector;
        Uvector= SeaStates.Uvector;
        p= SeaStates.p;
        seaType= SeaStates.seaType;
    else
        % DEFAULT SEA STATE 6
        Hvector= 3;
        wpVector= 2*pi/16;
        Uvector= 16;
        p= 1;
        seaType= 'Bret';
    end
    
else %Defaults
    geoMod= 'none';
    cylParam= [];
    TMDeParam= [];
    platformParam= [];
    heavePlate= [];
    frameParam= [];
    nonLinParam= [];
    Cd= 0;
    mSwitch= 'lin'; %'lin' or 'Nonlin'
    cSwitch= 'lin';
    domainSwitch= 'frequency'; % 'time';
    runType= 'iterate'; %'RAO'
    %fp= .0666;
    WellsParam= [0 0];
    useWpAB= 0; %indicator to use constant A and B matrices
end

%evalin('caller','mfilename')

%% PLOT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotPlatformTowerMotion= 0;
plotPdf= 0;
plotOWC=0; %Set to 1 to output OWC plots
barPlots= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TOWER MATERIAL PROPERTIES

% Used in calculate_Fatigue.m, which is called by simulateFWT_OWC_sphere_OmegaBeta_20160908
sigma_ult= 2260e6; %415e6; %5420e6; %[Pa]
m= 5; %4; %5.5; 
% n_Life= 0; %initially set number of fatigue cycles to 0. calculate_Fatigue.m adds to n_Life
materialParam= [sigma_ult, m];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
%% ITERATE THROUGH WIND/SEA STATES
Stress= zeros(length(p), 1); %in [Pa]
StressMean= zeros(length(p), 1); %in [Pa]
X1= zeros(length(p), 1);
X2= zeros(length(p), 1);
X3= zeros(length(p), 1);
X4= zeros(length(p), 1);
X5= zeros(length(p), 1);
X6= zeros(length(p), 1);
X1N1= zeros(length(p), 1); %surge of nacelle relative to platform
X2N1= zeros(length(p), 1); %sway of nacelle relative to platform
X1N2= zeros(length(p), 1); %surge of nacelle relative to platform
X2N2= zeros(length(p), 1); %sway of nacelle relative to platform
XTMD= zeros(length(p), 1);
XTMDe= zeros(length(p), 1);
ZTMDe= zeros(length(p), 1);
XN= zeros(length(p), 1); %sway of nacelle relative to platform
V1= zeros(length(p), 1); %platform velocities
V2= zeros(length(p), 1);
V3= zeros(length(p), 1);
V4= zeros(length(p), 1);
V5= zeros(length(p), 1);
V6= zeros(length(p), 1);
V1N1= zeros(length(p), 1); %surge velocity of nacelle due to bending
V2N1= zeros(length(p), 1); %sway velocity of nacelle due to bending
V1N2= zeros(length(p), 1); %surge velocity of nacelle due to bending
V2N2= zeros(length(p), 1); %sway velocity of nacelle due to bending
arrayPower= zeros(length(p), 1);
arrayPower2= zeros(length(p), 1);
arrayZ= zeros(length(p), 1);
tubeZ= zeros(length(p), 1);
dVector= zeros(length(p), 1);
D_towerBase_surge= zeros(length(p), 1);
n_stateVector= zeros(length(p), 1);
PMass= zeros(length(p), 1);
% XssOut= zeros(16, 1);
XssVector= zeros(40,length(p));
VssVector= zeros(40,length(p));

windPowers= zeros(length(p), 1);

for i= 1:length(p)
    i
    
    [D_state, n_state, Stress_Tower, stress_mean, Xss, Vss, ap, ap2, V33_WEC_wrt_Platform_i, XNacelle, windPower]=...    
    simulateFWT_WEC_20190519(runType, seaType, wpVector(i), Hvector(i),...
    Uvector(i), Cd, mSwitch, cSwitch, domainSwitch, p(i), geoMod, materialParam,...
    cylParam, TMDeParam, platformParam, heavePlate, frameParam, WellsParam, useWpAB); 

    Stress(i)= Stress_Tower; %RMS stress amplitude [Pa]
    StressMean(i)= stress_mean; %[Pa]
    X1(i)= Xss(1); %RMSs
    X2(i)= Xss(2);
    X3(i)= Xss(3);
    X4(i)= Xss(4);
    X5(i)= Xss(5);
    X6(i)= Xss(6);
    X1N1(i)= Xss(7);
    X2N1(i)= Xss(8);
    X1N2(i)= Xss(9);
    X2N2(i)= Xss(10);
    XTMD(i)= Xss(13+ sum(strcmp(geoMod, 'TMD')));
    XTMDe(i)= Xss(13+ sum(strcmp(geoMod, 'TMD')) +sum(strcmp(geoMod, 'TMDe')) );
    ZTMDe(i)= Xss(13+ sum(strcmp(geoMod, 'TMD')) +2*sum(strcmp(geoMod, 'TMDe')) );
    XN(i)= XNacelle;
    V1(i)= Vss(1);
    V2(i)= Vss(2);
    V3(i)= Vss(3);
    V4(i)= Vss(4);
    V5(i)= Vss(5);
    V6(i)= Vss(6);
    V1N1(i)= Vss(7);
    V2N1(i)= Vss(8);
    V1N2(i)= Vss(9);
    V2N2(i)= Vss(10);
%     arrayZ(i)= Zowc; %RMS displacement of OWC with maximum amplitude
    arrayPower(i)= ap; %total array power [W]
    arrayPower2(i)= ap2; %total array power [W]
%     tubeZ(i)= tz; %RMS displacement of tube with maximum amplitude
%     dVector(i)= dOut;
    V33_WEC_wrt_Platform(i) = V33_WEC_wrt_Platform_i;

    D_towerBase_surge(i)= D_state;
    n_stateVector(i)= n_state;
    XssOut(:,i)= Xss;
    windPowers(i)= windPower;
end

%% DAMAGE, DEL, AND FATIGUE STRESS

D_tot= sum(D_towerBase_surge);%.*p);
    TlifeYr= 20; %Design lifetime in seconds
    D_oneYear= D_tot./TlifeYr;
sigma_m_Life= sum(StressMean.*p'); %[Pa] %average mean stress during lifetime
n_Life= sum(n_stateVector);
DEL_stress= 2*( D_tot/(2*n_Life) )^(1/m)*(sigma_ult-abs(sigma_m_Life)); %this DEL is the pk-pk stess amplitude [Pa]

IBs= 2.876; %[m^4]
roBs= 3.25; %[m]
Kc= 1.6; %Concentration factor beyond normal beam bending. See notebook 13 pages 39, 147
% TwrBsMMxy_kNm= abs(Stress).*IBs./roBs./10^3./Kc;    %tower base bending moment based on beam. 0-pk MONOCHROMATIC ONLY [kNm]

DEL= DEL_stress.*IBs./roBs./10^3./Kc;    %tower base bending moment based on beam stress pk-pk [kNm]
maxStress= max(Stress); %maximum stress RMS that occurs during lifetime. NOTE: stats don't tell us maximum stress. [Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEAN MOTIONS AND STRESS 
% %nacelle velocity !!Correct these before use!!
% Vn_surge_total= abs(V1+L_SWL_to_Nacelle.*V5+V1N1); %platform surge + length*platform pitch + tower bending are complex amplitudes
% Xn_surge_total=  abs(X1+L_SWL_to_Nacelle.*X5+X1N1);
Pmat= repmat(p,length(XssOut(:,1)), 1);
Xssmean= sum(XssOut.*Pmat, 2);
Xssmax= max(XssOut,[], 2); %max of each row
% X5mean= sum(X5.*p');
% XTMDmean= sum(XTMD.*p');b
% XTMDemean= sum(XTMDe.*p');
% ZTMDemean= sum(ZTMDe.*p');
% X_surge_SWL= X1;
% 
% %velocity of platform at SWL (important for maintenance workers?)
% Vplatform_surge_total= abs(V1); %platform surge
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OWC ARRAY POWER AND COST %%%%%%%%%%%%%%%%%%%%%%%%%%%
arrayPower= arrayPower./10^3; %[KW]
arrayPower2= arrayPower2./10^3; %[KW]
meanPower= sum(arrayPower.*p'); %Power from first array [kW]
meanPower2= sum(arrayPower2.*p'); %Power from second array[kW]
% MeanPower = mean_V33_WEC_wrt_Platform_RMS * dJz * eta
% (eta= 0.8)
mean_V33_WEC_wrt_Platform_RMS= sum(V33_WEC_wrt_Platform'.*p'); % Realtive velocity between FWT and WEC [RMS]
% 
maxPower_KW= max(arrayPower); %maximum power from sea states [kW]
CF= meanPower/maxPower_KW;
% % CF(CF<0.3)= 0.3;
% maxPower_KW= meanPower/0.3;
if CF<0.3
    ii= length(arrayPower);
    [arrayPowerSort,indSort] = sort(arrayPower);
    pSort= p(indSort);
    while CF<0.3     
        arrayPowerSort(ii:end)= arrayPowerSort(ii-1);
        maxPower_KW= max(arrayPowerSort); %maximum power from sea states. Should be ii-1. [kW]
        meanPower= sum(arrayPowerSort.*pSort'); %Power from first array [kW]
        CF= meanPower/maxPower_KW;
        ii= ii-1;
    end
end


% Cap= max(arrayPower)./(10^6); %array capacity is maximum power over year. Already accounts for N tubes [MW]
% hinges= 0;end

% arrayCost= Fixed_Shallow_Cost(Cap, Dw, T_col, N, hinges, R_col);
% capFactor= meanPower./Cap;
% LCOEout= OWC_Fixed_Shallow_LCOE (Cap, capFactor, arrayCost); %[$/kWh]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WIND POWER %%%%%%%%%%%%%%%%%%%%%%%%%%%

windPowerOut= sum(windPowers.*p');
maxWindPower= max(windPowers);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
















