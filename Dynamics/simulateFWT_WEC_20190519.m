function [D_state, n_state, Stress_Tower, sigmaMean, X_RMS, V_RMS,...
    array_Power, array_Power2, V33_WEC_wrt_Platform, sigma_Nacelle, windPower, Hinge]=...
    simulateFWT_WEC_20190519(runType, seaType, w_p,...
    Hs, U, Cd, mSwitch, cSwitch, domainSwitch, p, geoMod, materialParam, cylParam,...
    TMDeParam, platformParam, heavePlate, frameParam, WellsParam, useWpAB)
% Creates dynamics matrix and computes responses in a sea state.
% For external TMDe: It is rigidly connected to FWT except for heave direction.
%                    Array is symmetric (ignore off-diagonal terms?).
%                    WEC properties are wrt SWL.

%% print RAO output parameters
plotRAO= 0;

%% CODE RUN CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(runType, 'RAO')
    RAOrun= 1; %Do not calculate fatigue. Do define the RAO output parameters
elseif strcmp(runType, 'iterate')
    RAOrun= 0;
elseif strcmp(runType, 'compute_hinges')
    RAOrun= 0;
end

%% WIND CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flags determine if use Linearized or stat Lin values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROTOR CONTROL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(domainSwitch, 'time')
    [angle_Kp, KpTable, angle_Ki, KiTable]= Jonkman_bladeGain; %angles in Deg.
    %NOTE: Torque-rate limit of 15,000 Nm/s
    [omegaT, Tgen, gradTgen]= genTorqueController;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOMETRY MODIFICATION INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values

%External TMD- may take the form of a cylinder, sphere, internal/external tuned mass damper.
TMDe= 0; %indicator of if 1 or more external(or internal) TMD's exist
numTMDe= 0;
zTMD= 0;

%Nonlinear TMDe
kjxn= 0;
kjyn= 0;
kzn= 0;
djxn= 0;
djyn= 0;
djzn= 0;

%Modify Platform- these parameters aren't used if modifyPlatform= 0
Ltop= 14;
LtopSub= 4;
Lmid= 8;
Lbot= 108;
Dtop= 6.5;
Dbot= 9.4;
zCM= -89.9155; %depth of platform CM in SWL coordinate system [m]
mPlatform= 7.466330e6; %platform mass
IplatformCMx= 4229230000; %[Kgm^2]
C33matplat= 332941; %[N/m]
C55matplat= 6.2795e+09 -4999180000; %[Nm/rad]
% mPlatformSteel= 1.7e6;
% mPlatformConc= 5.7e6;
% mSphereSteel= 0;
% tst= .0522;
IplatformCMzz= 164230000; %yaw inertia
modifyPlatform= 0;
platName= 'OC3';

if sum(strcmp(geoMod, 'Cylinder'))
    rCyl= cylParam(1);
    lCyl= cylParam(2);
    zCyl= cylParam(3);
    mCyl= cylParam(4); %NOT NECESSARILY neutrally buoyant [Kg]
    VCylSub= cylParam(5);
    rCyl1=9.4/2; %inner. OC3-Hywind
    rCyl2= rCyl;
    %     VCyl= pi*(rCyl-rCyl1)^2*lCyl;
    Cyl= 1; %indicator
else
    zCyl= 0;
    Cyl= 0;
end

if sum(strcmp(geoMod, 'heavePlate'))
    usePlate= 1;
    blah= 1;
else
    usePlate= 0;
end

if sum(strcmp(geoMod, 'TMDe'))
    
    %Each column is a different member of array
    
    %Geometry
    mTMDe= TMDeParam.mTMDe;
    I55TMDe= TMDeParam.I55TMDe;
    
    A11TMDe= TMDeParam.A11TMDe;
    A15TMDe= TMDeParam.A15TMDe;
    A51TMDe= TMDeParam.A51TMDe;
    A55TMDe= TMDeParam.A55TMDe;
    A33TMDe= TMDeParam.A33TMDe;
    
    C33TMDe= TMDeParam.C33TMDe;
    C55TMDe= TMDeParam.C55TMDe;
    
    %Wells parameters
    etaW= 0.6; %TMDeParam{5};
    dWellsx= TMDeParam.dWellsx;
    
    %Pressure EOM
    kWells= TMDeParam.kWells;
    VChamber= TMDeParam.VChamber;
    
    %Water column EOM
    ICol= TMDeParam.ICol; %includes added mass
    %C33Col= TMDeParam{10}; eom is normalized by area, so just rho*g used
    ACol= TMDeParam.ACol;
    
    %Parameters that may differ for each component of array
    %Each of these may be a horizontal vector
    dJx= TMDeParam.dJx;
    dJy= TMDeParam.dJy;
    dJz= TMDeParam.dJz;
    KJx= TMDeParam.KJx;
    KJy= TMDeParam.KJy;
    KJz= TMDeParam.KJz;
    etaJ= TMDeParam.etaJ;%.85; %.82;
    
%     K1J= TMDeParam.K1In;
%     K2J= TMDeParam.K2In;
%     K3J= TMDeParam.K3In;
%     d1J= TMDeParam.d1In;
%     d2J= TMDeParam.d2In;
%     d3J= TMDeParam.d3In;
    
    %Location of junction point on FWT
    x_junction_FWT= TMDeParam.Lfx;
    y_junction_FWT= TMDeParam.Lfy;
    z_junction_FWT= TMDeParam.Lfz;
    
    %Location of link end on WEC-side. May be in nacelle
    x_junction_WEC= TMDeParam.Lwx;
    y_junction_WEC= TMDeParam.Lwy;
    z_junction_WEC= TMDeParam.Lwz;
    zTMD= z_junction_WEC(1); %IF INTERNAL TURNED MASS DAMPER IN NACELLE- FIRST ELEMENT
    
    %Link vector. Each column is a different array component
    x_FWT_to_WEC_vec= x_junction_WEC - x_junction_FWT;
    y_FWT_to_WEC_vec= y_junction_WEC - y_junction_FWT;
    z_FWT_to_WEC_vec= z_junction_WEC - z_junction_FWT;
    
    % TMDe lateral and vertical forcing- depends on depth:
    
    % 2. For Surge Force: Input vectors of zS, xS, dVS, dAS
    zSVector= TMDeParam.zSVector;
    xSVector= TMDeParam.xSVector;
    dVSVector= TMDeParam.dVSVector;
    dASVector= TMDeParam.dASVector;
    
    %3. For Heave Force: Input vector of  zH, xH, dVH, dAH <-- for strip theory OR Sphere
    zHVector= TMDeParam.zHVector;
    xHVector= TMDeParam.xHVector;
    dVHVector= TMDeParam.dVHVector;
    dAHVector = TMDeParam.dAHVector;
    
    xWellsForceSVector= TMDeParam.xWellsForceSVector; %for water forcing hitting horizontal Wells
    zWellsForceSVector= TMDeParam.zWellsForceSVector; %for water forcing hitting Wells
    xWellsForceHVector= TMDeParam.xWellsForceHVector; %for water forcing hitting Wells
    zWellsForceHVector= TMDeParam.zWellsForceHVector; %For water forcing hitting bottom of column
    
    %     controlKJxHeave= TMDeParam{36};
    %     controlKJySurge= TMDeParam{37};
    %     controlKJyHeave= TMDeParam{38};
    %     controlKJzSurge= TMDeParam{39};
    
    
    kzn= TMDeParam.kzn; %WEC hydrostatic stiffness
    rigidPitch= TMDeParam.rigidPitch;
    %rigidTranslate= TMDeParam.rigidTranslate;
        
    dWellsz= TMDeParam.dWellsz;
    
    controlKxx= TMDeParam.controlK11;
    controlKzz= TMDeParam.controlK33;
    
    TMDe= 1; %indicator of if external TMD exists
    
    numTMDe= length(x_junction_WEC);
    
    %%
    if controlKxx==1 %
        k11TMDe= (mTMDe+A11TMDe).*w_p^2;% - C33TMDe;
        k11TMDe= max(k11TMDe, 0); %Restrict to greater than 0
        K1J= k11TMDe;
    end
    if controlKzz==1
        k33TMDe= (mTMDe+A33TMDe).*w_p^2;% - C33TMDe;
        k33TMDe= max(k33TMDe, 0); %Restrict to greater than 0
        K3J= k33TMDe;
    end
    
end

if sum(strcmp(geoMod, 'Platform'))
    platformParam1= platformParam;%{1};
    zCM= platformParam1(1); % used for mass terms
    LtopSub= platformParam1(2);
    mPlatform= platformParam1(3);
    IplatformCMx= platformParam1(4);
    C55matplat= platformParam1(5);
    C33matplat= platformParam1(6);
    %     mPlatformSteel= platformParam1(7);
    %     mPlatformConc= platformParam1(8);
    %     mSphereSteel= platformParam1(9);
    Ltop= platformParam1(10);
    Lmid= platformParam1(11);
    Lbot= platformParam1(12);
    Dtop= platformParam1(13);
    Dbot= platformParam1(14);
    %     tst= platformParam1(15);
    modifyPlatform= 1;
    platName= platformParam1(15); %OC4 or blank
end
if sum(strcmp(geoMod, 'Tower'))
    
end

if sum(strcmp(geoMod, 'frame'))
    frameM= frameParam{1};
    frameC= frameParam{2};
else
    frameM= zeros(6,6);
    frameC= zeros(6,6);
end

dWells= WellsParam(1);
zWells= WellsParam(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEA CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(platName, 'OC4')
    H= 200;
else
    H= 320; %water depth [m]
end

[Su, wVector, kVector, VgVector, dw]= declareSeaConditions(seaType, w_p, Hs, H);
% Allowable inputs for seaType: 'Bret', 'Jon', 'white', 'mono'

g= 9.81; %acceleration due to gravity [m/s^2]
rho= 1025; %sea water density [Kg/m^3]

gamma= 1.4; %air specific heat ratio
Patm= 101325; %[Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLATFORM GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coordinate system origin is at z= SWL, aligned with platform/tower longitudinal axis.

% Used for viscous damping and mooring laod at fairleads
if strcmp(platName, 'OC4')
    %Don't change viscous parameters for now
    zDraft= -(Ltop+Lmid+Lbot);%[m]
    zzP= linspace(zDraft, 0, 10); %submerged coordinates
    dzzP= zzP(2)-zzP(1);
    D= Dbot.*ones(1,length(zzP));
    D(zzP>(-Ltop-Lmid))= (Dtop-Dbot)/Lmid.*(zzP(zzP>(-Ltop-Lmid))+Ltop)+Dtop;
    D(zzP>-Ltop)= Dtop; %Diameter for viscous damping
    L_platform_raised= Ltop-LtopSub; %distance that platform extends above SWL [m]
    
    Zfl= -14; %z coordinate of mooring fairleads relative to still water level [m]
elseif strcmp(platName, 'OC3') %Assume spar platform
    zDraft= -(Ltop+Lmid+Lbot);%[m]
    zzP= linspace(zDraft, 0, 10); %submerged coordinates
    dzzP= zzP(2)-zzP(1);
    D= Dbot.*ones(1,length(zzP));
    D(zzP>(-Ltop-Lmid))= (Dtop-Dbot)/Lmid.*(zzP(zzP>(-Ltop-Lmid))+Ltop)+Dtop;
    D(zzP>-Ltop)= Dtop; %Diameter for viscous damping
    L_platform_raised= Ltop-LtopSub; % distance that platform extends above SWL [m]
    % used for eigen calculation
    
    Zfl= -70; %z coordinate of mooring fairleads relative to still water level [m]
else %if strcmp(platName, 'tapered')
    % 5/25/2019
    zDraft= -Lbot;%[m]
    zzP= linspace(zDraft, 0, 10); %submerged coordinates
    dzzP= zzP(2)-zzP(1);

    Dslope= ( Dbot - Dtop ) / -abs(Lbot);
    D= Dtop + Dslope.*zzP; %Diameter for viscous damping
    L_platform_raised= Ltop-LtopSub; % distance that platform extends above SWL [m]
    % used for eigen calculation
    
    Zfl= -45;
end

%% ROTOR GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Drotor= 125.88; %rotor diameter from Jonkman [m]
% Rrotor= Drotor/2;
Rrotor= 62.94*cosd(-2.5);
rho_a= 1.225; %[Kg/m^3]
S= pi*Drotor^2/4;

Ihub= 115.926e3; %[Kg m^2]
Irotor= 3*11.776047e6; %[Kg m^2]
Ngear= 97;
Igen= 534.116; %[Kg m^2]
Idrivetrain= Ihub+Irotor+Ngear^2*Igen;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET WAMIT COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if modifyPlatform
    if strcmp(platName, 'OC4')
        [wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable, A11Table, A15Table, A22Table, A24Table, A33Table, A44Table, A55Table, A66Table, B11Table, B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table]...
            = getWAMIToutputsOC4; %OC4
    elseif strcmp(platName, 'OC3')
        [wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable, A11Table, A15Table, A22Table, A24Table, A33Table, A44Table, A55Table, A66Table, B11Table, B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table]...
        = getWAMIToutputs; %OC3
    %     [A11INF, A15INF, A22INF, A24INF, A33INF, A42INF, A44INF, A51INF, A55INF, A66INF]= getWAMIT_infinite_coeffs;
    else
        [wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable, A11Table, A15Table, A22Table, A24Table, A33Table, A44Table, A55Table, A66Table, B11Table, B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table]...
            = predictplatformLateralForces(H, zCM, zDraft, zzP, dzzP, D);
    end
else % OC3
     [wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable, A11Table, A15Table, A22Table, A24Table, A33Table, A44Table, A55Table, A66Table, B11Table, B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table]...
        = getWAMIToutputs; %OC3
    %     [A11INF, A15INF, A22INF, A24INF, A33INF, A42INF, A44INF, A51INF, A55INF, A66INF]= getWAMIT_infinite_coeffs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp(geoMod, 'Tower'))
    mTower
    mNacelleHub
    mRotor
    L_SWL_to_Nacelle
    IzzHubNacelle % I66
    IzzTower
    
    mMooring= 3*(77.1*902.2);
    
    [K77, M77, K99, M99, M71, M75, M91, M95, K79, M79,...
        m11_tower, m55_tower, m15_tower, m51_tower,...
        B77T, B17T, B71T, B75T, B11T, B55T, B15T, B99T,...
        u1TMD, u2TMD]=...
        normalize_tapered_Eigen_DUMMY_20190213(zTMD, L_platform_raised);
    
else
    %% TOWER GEOMETRY A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E= 210e9; %steel elastic modulus [Pa]
    rho_st= 8500; %structure density [Kg/m^3]
    mTower= 249718; %347460; %[Kg]
    mNacelleHub= 56780+240000; %hub + nacelle [Kg]
    mRotor= 3*17740; %rotor mass
    mMooring= 3*(77.1*902.2);
    % mEnd= 3.5e5;
    %% TOWER GEOMETRY B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L_SWL_to_Nacelle= 87.6;
    
    %Izz
    IzzHubNacelle= 2607890+115926; %nacelle+hub [Kg m^2]
    
    z= linspace(0, 77.6, 300); %tower coordinates [m]
    dz= z(2)-z(1);
    ro=  3.25 - .0169.*z; %outer radius [m]
    ri= 3.223- .0168.*z; %inner radius [m]
    dm= rho_st.*dz.*pi.*(ro.^2-ri.^2);
    dIz= .5.*dm.*(ri.^2 + ro.^2);
    IzzTower= sum(dIz);
    
    [K77, M77, K99, M99, M71, M75, M91, M95, K79, M79,...
        m11_tower, m55_tower, m15_tower, m51_tower,...
        B77T, B17T, B71T, B75T, B11T, B55T, B15T, B99T,...
        u1TMD, u2TMD]=...
        normalizeEigenShape_Jonkman_really_20160803(zTMD, L_platform_raised);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%
numDOF= 10+3+3*numTMDe; %6 platform + 4 tower + Omega(11) + qDotw(12) + deltaBeta(13) + zTMDe(14) + zCol(15) + Pc(16) 12/31/2018
% firstOWCDOF= 10+3+1;
%% 1/w COEFF MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL= zeros(numDOF, numDOF);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MASS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M= zeros(numDOF, numDOF);

% %Platform
M(1,1)= (mPlatform + m11_tower); %mTower + mPlatform+mNacelleHub; %surge
M(2,2)= mPlatform+ m11_tower; %mTower + mPlatform+mNacelleHub; %sway
M(3,3)= mPlatform + mTower + mNacelleHub+ mRotor +mMooring; %heave
M(4,4)= IplatformCMx + mPlatform*zCM^2 + m55_tower; % [kgm^2]
M(5,5)= (IplatformCMx + mPlatform*zCM^2 + m55_tower); %[kgm^2]
M(6,6)= IplatformCMzz+ IzzTower+IzzHubNacelle; %[Kgm^2]

% M(1,5)= mPlatform*zCM + 1.5*m15_tower;
% M(5,1)= mPlatform*zCM + 1.5*m51_tower;
M(1,5)= (mPlatform*zCM + m15_tower);
M(5,1)= (mPlatform*zCM + m51_tower);

M(2,4)= -mPlatform*zCM - m15_tower;
M(4,2)= -mPlatform*zCM - m51_tower;


% %EXPERIMENT PLATFORM
% M(4,4)= 6.9e10; %IplatformCMx + mPlatform*zCM^2 + m55_tower; % [kgm^2]
% M(5,5)= 6.9e10; %IplatformCMx + mPlatform*zCM^2 + m55_tower; %[kgm^2]
% M(1,5)= 2.1417e5; %mPlatform*zCM + m15_tower;
% M(5,1)= 2.1417e5; %mPlatform*zCM + m51_tower;
% M(2,4)= -2.1417e5; %-mPlatform*zCM - m15_tower;
% M(4,2)= -2.1417e5; %-mPlatform*zCM - m51_tower;


M(1,7)= M71;
M(5,7)= M75;
M(1,9)= M91;
M(5,9)= M95;

%Tower bending equations
M(7,7)= M77; %mode 1, fore-aft
M(7,1)= M71;
M(7,5)=  M75;
M(8,8)= M77;  %Mode 1, side-side
M(8,2)= M71;
M(2,8)= M71;
M(8,4)= -M75;
M(4,8)= -M75;
M(7,9)= M79;
M(9,7)= M79;

M(9,9)= M99; %mode 2, fore-aft
M(9,1)= M91;
M(9,5)=  M95;
M(10,10)= M99; %mode 2, side-side
M(10,2)= M91;
M(2,10)= M91;
M(10,4)= -M95;

%Blade pitch equation
% M(13,13) = 1; Use C(13,13)= 1.

%CYLINDER
if Cyl
    Izcyl= mCyl*(rCyl1^2+rCyl2^2); %about z
    IxCyl= mCyl/12*(3*(rCyl1^2+rCyl2^2) + lCyl^2);
    M(1:6, 1:6)= M(1:6, 1:6) + [mCyl            0       0           0       mCyl*zCyl        0 %zSphere is negative. See Notebook 15 p. 157
        0            mCyl       0     -mCyl*zCyl        0            0
        0               0       mCyl        0           0            0
        0        -mCyl*zCyl      0  IxCyl+mCyl*zCyl^2   0            0
        mCyl*zCyl       0        0          0    IxCyl+mCyl*zCyl^2   0
        0               0        0          0           0           Izcyl];
end

% HEAVE PLATE
if usePlate
    
    M(1,1)= M(1,1) + heavePlate.mTMDe;
    M(1,5)= M(1,5) + heavePlate.mTMDe * heavePlate.Lwz;
    M(3,3)= M(3,3) + heavePlate.mTMDe;
    M(5,1)= M(5,1) + heavePlate.mTMDe * heavePlate.Lwz;
    M(5,5)= M(5,5) + heavePlate.I55TMDe + heavePlate.mTMDe * heavePlate.Lwz^2;
    
end


%External TMD
if TMDe
    % r is distance from y axis: value of x.
    r= x_junction_WEC;%( Lwx.^2 + Lwy.^2 ).^.5;
    for ii= 1:numTMDe
        %         xwInd= 14+5*ii;
        %         ywInd= 15+5*ii;
        %         M(xwInd,xwInd)= M(xwInd,xwInd) + mTMDe;
        %         M(ywInd,ywInd)= M(ywInd,ywInd) + mTMDe;
        
        % See thesis PDF page 58. zCol is relative motion
        
        % add rigid properties to FWT platform:
        % only I55 accounts for water column.
        % A11, A15 account for water column.
        I55_new= I55TMDe + mTMDe*r(ii)^2;
        m_eff= mTMDe;
        
        % ignore I66
        M(1:6,1:6)= M(1:6,1:6) + ...
            [m_eff    0      0       0       0       0
            0     m_eff   0       0       0       0
            0       0     0       0       0       0
            0       0     0    I55_new    0       0
            0       0     0       0    I55_new    0
            0       0     0       0       0       0];
        
        % WEC's free properties:
        zwInd= 11+3*ii;
        zColInd= 12+3*ii;
        
        M(zwInd,zwInd)= M(zwInd,zwInd) + mTMDe;
        mCol= ICol/ACol;
        mCol(isnan(mCol))= 0;
        M(zColInd,zColInd)= M(zColInd,zColInd) + mCol; % Account for column added mass here
        M(zColInd,zwInd)= M(zColInd,zwInd) + mCol; % zCol is relative to WEC heave
    end
end

%FRAME
M(1:6,1:6)= M(1:6,1:6) + frameM;

if exist('Hinge') && isStruct(Hinge)
    for ii= 1:numTMDe
        M(1,1)= M(1,1) + Hinge.mass; % Hinge contains: mass, r_out, r_in, L
        M(2,2)= M(2,2) + Hinge.mass;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONSTANT ADDED MASS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A= zeros(numDOF,numDOF); %wind turbine values added during iterations

% Cylinder
if Cyl
    A11Cyl= rho*VCylSub/2;
    % A11Cyl= A11Cyl + rho*VCylSub/2; %0 frequency limit [Kg]
    A22Cyl= A11Cyl;
    A33Cyl= A11Cyl;
    
    A15Cyl= A11Cyl*zCyl; %zSphere is negative. See Notebook p. 157. Code below duplicates for A51
    A24Cyl= -A11Cyl*zCyl; %assume rotational symmetry
    A44Cyl= A11Cyl*zCyl^2;
    A55Cyl= A44Cyl;
    A66Cyl= 0;
else
    A11Cyl= 0;
    A22Cyl= A11Cyl;
    A33Cyl= A11Cyl;
    
    A15Cyl= 0; %zSphere is negative. See Notebook p. 157. Code below duplicates for A51
    A24Cyl= 0; %assume rotational symmetry
    A44Cyl= 0;
    A55Cyl= 0;
    A66Cyl= 0;
end


% HEAVE PLATE
if usePlate
    
    A(1,1)= A(1,1) + heavePlate.A11TMDe;
    A(1,5)= A(1,5) + heavePlate.A11TMDe * heavePlate.Lwz;
    A(3,3)= A(3,3) + heavePlate.A33TMDe;
    A(5,1)= A(5,1) + heavePlate.A11TMDe * heavePlate.Lwz;
    A(5,5)= A(5,5) + heavePlate.A55TMDe + heavePlate.A11TMDe * heavePlate.Lwz^2;
    
end

% external tuned mass damper
if TMDe
    for ii=1:numTMDe
        A24TMDe= A15TMDe;
        % off-diagonal term will be cancelled out
        % due to array symmetry:
        % A15_eff= A15TMDe + Lwx(ii)*A11TMDe
        % sum for dz [-L to 0]:
        A55_new= A55TMDe(ii) + A11TMDe(ii)*x_junction_WEC(ii)^2;
        
        % add WEC rigid properties to FWT platform:
        % I55 accounts for water column.
        % A11, A15 account for water column.
        A(1:6,1:6)= A(1:6,1:6) + ...
            [A11TMDe(ii)         0    0       0       A15TMDe(ii) 0
            0        A11TMDe(ii)   0    A24TMDe(ii)    0       0
            0            0         0       0           0       0
            0        A24TMDe(ii)   0    A55_new        0       0
            A51TMDe(ii)     0         0       0        A55_new    0
            0            0         0       0           0       0];
        
        % WEC's free properties:
        zwInd= 11+3*ii;
        zColInd= 12+3*ii;
        A(zwInd, zwInd)= A33TMDe(ii);
        
        % water column added mass is already accounted for in mCol
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTANT DAMPING MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B= zeros(numDOF, numDOF); %Correction for adjusting B while iterating.

% Linear wind damping defined after determine mean Omega. Added to B matrix

%Tower damping: B77T B17T B71T B57T B75T B11T B55T B15T
B(7,7)= B77T;
B(1,7)= B17T;
B(7,1)= B71T;
B(5,7)= B75T;
B(7,5)= B75T;
B(1,1)= B11T;
B(5,5)= B55T;
B(1,5)= B15T;
B(5,1)= B15T;
% b_towerBend1= .01.*2; % 1% damping ratio
% b_towerBend2= 0; %.01*2*(K_towerBend2*m_TowerBend2)^.5; % 1% damping ratio

B(11,11)= Idrivetrain; %Omega0 eqn.

%qDotw equation
B(12,12)= 1;
B(12, 1)= 1;
B(12, 5)= L_SWL_to_Nacelle;

%Beta equation
% B(13,13)= 2*.02*30*2*pi; %206000*pi/180;

%WATER COLUMN ARRAY damping and other dampings populated within w-iteration
%Tube damping
% columnDampingOptimization; %outputs optimized damping coefficient of each column. Adds damping to B.

%External TMD- Rotational Dampers
% junction_damping;


% Wells turbine rigidly attached to spar
B(1,1)= B(1,1) + dWells;
B(1,5)= B(1,5) + dWells*zWells;
B(5,1)= B(5,1) + dWells*zWells;
B(5,5)= B(5,5) + dWells*zWells^2;

if TMDe
    for ii= 1:numTMDe
        % WEC's free properties:
        zwInd= 11+3*ii;
        zColInd= 12+3*ii;
        PInd= 13+3*ii;
        
        B(PInd, PInd)= kWells*VChamber/(gamma*Patm);
        B(PInd, zColInd)= -kWells*ACol;
        
        % Junction between FWT and WEC (1st order approximation)       
        B(3, 3)= dJz;
        B(3, zwInd)= -dJz;
        B(zwInd, 3)= -dJz;
        B(zwInd, zwInd)= dJz;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C= zeros(numDOF,numDOF);
C(3,3)= C33matplat;
C(4,4)= C55matplat;
C(5,5)= C55matplat;

%Hydrostatic

% %EXPERIMENT
% C(3,3)= 3.255e5; %3.35e5; % C33matplat;
% C(4,4)= 2.8065e8; %1.18e9; %C55matplat;
% C(5,5)= 2.8065e8; %1.18e9; %C55matplat;

%Mooring lines- linearized stiffness determined in DC section


%Tower stiffness
C(7,7)= K77;
C(8,8)= K77;
C(9,9)= K99;
C(10,10)= K99;
C(7,9)= K79;
C(9,7)= K79;

%qdotW equation
C(12,12)= 1;

%Beta equation
C(13,13)= 1;
%C(13,13)= 30*2*pi; %971350000*pi/180;

%CYLINDER
if Cyl
    if zCyl==0 %Is sphere at surface, modify
        C(3,3)= C(3,3) + rho*g*pi*(rCyl2^2-rCyl1^2);
        IxxCyl= pi/4*(rCyl2^4-rCyl1^4);
        zBarSub= zCyl; %NEED TO CHANGE IF EVER HAVE CYLINDER SWL
    else
        IxxCyl= 0;
        zBarSub= zCyl;
    end
    %Assume symmetric, by definition of cylinder
    %Even if completely submerged, changes rotational stiffness
    C(4,4)= C(4,4) + rho*g*IxxCyl + rho*VCylSub*zBarSub - mCyl*g*zCyl;
    C(5,5)= C(5,5) + rho*g*IxxCyl + rho*VCylSub*zBarSub - mCyl*g*zCyl;
end

% TUBE ARRAY
% tubeStiffness; %output matrix of tubes in OWC array- when rigidly attached


%EXTERNAL TMD
if TMDe
    
    for ii= 1:numTMDe
        
        % add rigid properties to FWT platform:
        % Array symmetry -> C35= 0.
        C55_new= C55TMDe + 0; % C11= 0. C11*Lwx= 0.
        C(1:6,1:6)= C(1:6,1:6) + ...
            [0     0    0       0         0      0
            0     0    0       0         0      0
            0     0    0       0         0      0
            0     0    0    C55_new      0      0
            0     0    0       0       C55_new  0
            0     0    0       0           0    0];
        
        % WEC's free properties:
        zwInd= 11+3*ii;
        zColInd= 12+3*ii;
        PInd= 13+3*ii;
        
        C(zwInd:PInd, zwInd:PInd)= [C33TMDe 0 -ACol
            rho*g rho*g 1
            0      0    1];
        
        % Junction between FWT and WEC (1st order approximation)       
        C(3, 3)= KJz;
        C(3, zwInd)= -KJz;
        C(zwInd, 3)= -KJz;
        C(zwInd, zwInd)= KJz;
        
        % These are only used to determine stress in the hinges.
        Lxl= x_FWT_to_WEC_vec(ii);
        Lyl= y_FWT_to_WEC_vec(ii);
        Lzl= z_FWT_to_WEC_vec;
        Llink= [Lxl Lyl Lzl];
        
        %Assume that junction point on FWT is in center of platform axis.
        Lxj= x_junction_FWT(ii);
        Lyj= y_junction_FWT(ii);
        Lzj= z_junction_FWT;
        
        
%         junction_springs; % Rotational Springs
%         translation_springs; % Rigid Links- or translation spring
        
    end
end

%FRAME
C(1:6,1:6)= C(1:6,1:6) + frameC;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FORCING VECTOR TABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tables from Jonkman

%Surge
X1= interp1(wForceTable, X1magTable, wVector); %NORMALIZED BY WAVE AMPLITUDE
X1phase= interp1(wForceTable, X1phaseTable, wVector);

%Heave
X3= interp1(wForceTable, X3magTable, wVector); %NORMALIZED BY WAVE AMPLITUDE
X3phase= interp1(wForceTable, X3phaseTable, wVector);

%Pitch
X5= interp1(wForceTable, X5magTable, wVector); %NORMALIZED BY WAVE AMPLITUDE
X5phase= interp1(wForceTable, X5phaseTable, wVector);

%Tables of complex excitation force amplitude at each frequency
F1= -X1.*exp(-1i.*X1phase); %NORMALIZED BY WAVE AMPLITUDE
F3= -X3.*exp(-1i.*X3phase); %NORMALIZED BY WAVE AMPLITUDE
F5= -X5.*exp(-1i.*X5phase); %NORMALIZED BY WAVE AMPLITUDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIAL VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
converge= 0;
it= 0; %number of iterations until converges

%covariance matrix. Row 1= Omega. Row 2= qDotw. Row 3: beta
Sigma= [.25 0 0; 0 .25 0; 0 0 .25]; %initial assumption for first iteration to run
sigmaQdot= zeros(length(zzP),1); %initial value for standard deviation of qDot
sigmaQ= 0;

%Initial values for nonlinear TMDe parameters
Sigma5JTMDe= 0;
% SigmaZrTMDe= 0;
% SigmaXrPTMDe= .1;
Sigma4JTMDe= .1;
SigmaWECheave= zeros(1,numTMDe);

%Borig defined after determine linearized wind damping
% Corig defined after define linear mooring forces below.

S_response= zeros(numDOF,length(wVector)); %square of response for each frequency of sea state
TF_Stress_surge= zeros(1,length(wVector)); %store tower bending amplitude
TF_terms= zeros(numDOF,length(wVector)); %transfer function values at each frequency
X_terms= zeros(numDOF,length(wVector)); %Force per amplitude on each DOF at each frequency
SigmaTemp= zeros(6,length(wVector));
F1bar_S= zeros(1,length(wVector));
M5bar_S= zeros(1,length(wVector));

%Simulation outputs
X_RMS=zeros(numDOF, 1);
X_RMS_old= zeros(numDOF, 1);
array_Power= 0; %initial value
% array_Power2= 0; %initial value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%% DC Values and mooring lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%Assume controls are about still-platform equilibrium

[windTableO, OmegaTable, windTableB, BetaTable]= JonkmanSS_Om_Beta; %steady-state values [OmegaTable is in rad/s. BetaTable is in deg]
OmegaMean= interp1(windTableO, OmegaTable, U); %[rad/s]
BetaMean1= interp1(windTableB, BetaTable, U);

[UbuoyTable, FthrustTable]= DCFthrustTable;
Fthrust= interp1(UbuoyTable, FthrustTable, U); %thrust in [N]
Mthrust= Fthrust*L_SWL_to_Nacelle;

if strcmp(platName, 'OC4')
    [q_table, F, dFdq, dq, c33LIN]= nonlinearMooringLinesOC4;
    % elseif strcmp(platName, 'OC3')
else %Assume OC3 mooring lines
    [q_table, F, dFdq, dq, c33LIN]= nonlinearMooringLines;
end

C55Mooringtemp= Zfl.^2.*interp1(q_table, -dFdq, 0);% + 109318000; %Negligible higher order effect- C55temp mooring compared to platform ~15%
ff= @(qq) [Fthrust + interp1(q_table, F, qq(1) ) %[Fthrust + interp1(q_table, F, qq(1) + Zfl.*qq(2)) %balance forces
    Mthrust + Zfl*interp1(q_table, F, qq(1) + Zfl.*qq(2)) - (C55matplat+C55Mooringtemp)*qq(2)]; %balance moments
%Mthrust + Zfl*interp1(q_table, F, qq(1) + Zfl.*qq(2)) -(C55matplat +Zfl.^2.*interp1(q_table, dFdq, qq(1) + Zfl.*qq(2))).*qq(2)];

[xDC, FVAL, EXITFLAG]= fsolve(ff, [0 0]); %qDC= [x1 x5 q]
x5mean= xDC(2); %[rad]
qDC= xDC(1) + Zfl*xDC(2); %mean displacement at mooring line height
C11m= interp1(q_table, -dFdq, xDC(1) + Zfl.*xDC(2));  %mooring line stiffness in neighborhood of DC value

% C11m= 0;%196875; %Check for experiment
if strcmp(mSwitch, 'lin') %Only applied to modes 1 and 5 for now
    %C11m= 41180; %for OC# linearized about 0 displacement
    %c15m= C11m*Zfl;
    %c55m= C11m*Zfl^2;
    C(1,1)= C(1,1) + C11m; %41180;         %June 20, 2016: Replace with nonlinear stiffness. 4.118e4
    C(1,5)= C(1,5) + C11m*Zfl; %-2882600; %Nov 29, 2016- use C11*L        -2821000;        %June 20, 2016: Replace with nonlinear stiffness
    C(5,1)= C(5,1) + C11m*Zfl; %-2882600; %Nov 29, 2016- use C11*L        % June 20, 2016: Replace with nonlinear stiffness
    C(5,5)= C(5,5) + C11m*Zfl^2; %+ 201782000; %3.11100000e8;  %Nov 29, 2016- use C11*L     % June 20, 2016: Replace with nonlinear stiffness
end
C(3,3)= C(3,3) + c33LIN; %11940;

%FOR OC3 for now- these modes uncoupled from 1,3,5
C(2,2)= C(2,2) + C11m; %+ 41180;
C(2,4)= C(2,4) - C11m*Zfl; %+2882600; %+ 2821000;
C(4,2)= C(4,2) - C11m*Zfl; %+2882600; %+ 2816000;
C(4,4)= C(4,4) + C11m*Zfl^2; %201782000; %311100000;
C(6,6)= C(6,6) + 11560000 + 9.834e7; %1.15e7


BetaMean= BetaMean1;%+ x5mean*180/pi; %For now, before add in new DOF [deg]
[qW_Matrix, Omega_Matrix, Ct_Matrix, Cq_Matrix, dCt_dq_Matrix, dCq_dq_Matrix, dCt_dOm_Matrix, dCq_dOm_Matrix, Cp_Matrix]= myThrust_versus_Urel_Om_Matrix_JK(BetaMean);

if strcmp(cSwitch, 'Nonlin')
    [qW_Matrix3, Omega_Matrix3, Beta_Matrix3, Ct_Matrix3, Cq_Matrix3, dCt_dq_Matrix3, dCq_dq_Matrix3, dCt_dOm_Matrix3, dCq_dOm_Matrix3, dCt_dB_Matrix3, dCq_dB_Matrix3, Cp_Matrix3]= myThrust_versus_Urel_Om_Matrix3_JK; %coefficients as 3D matrices (qW, Om, Beta)
    % %Each    row (x-value) is different Omega.
    % %Each column (y-value) is different qW.
    % %Each   page (z-value) is different beta.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LINEAR WIND DAMPING
if ~strcmp(cSwitch, 'Nonlin') %if linear
    %% LINEARIZED WIND DAMPING
    % [Vtable, CtTable]= fThrust_OC3_Baseline;
    % Ct= interp1(Vtable, CtTable, U);
    Cp= interp2(qW_Matrix, Omega_Matrix, Cp_Matrix, U, OmegaMean);
    Ct= interp2(qW_Matrix, Omega_Matrix, Ct_Matrix, U, OmegaMean);
    
    % linearization. Matches statistical linearization when U>>motion!
    b11w_Lin= rho_a*Ct*S*U;
    b15w_Lin= rho_a*Ct*S*U*L_SWL_to_Nacelle;
    b51w_Lin= b15w_Lin;
    b55w_Lin= rho_a*Ct*S*U*L_SWL_to_Nacelle^2;
    
    % %% Rotor power
    % % qWDot= X_RMS(12);
    %
    % windPower= .5*rho_a*S*Cp.*U.^3./10^6; %[MW]. Power= .5rhoA*S*E{Cp*qwDot^3}
    windPower= .5*rho_a*S*Cp.*U.^3./10^6; %[MW]. Power= .5rhoA*S*E{Cp*qwDot^3}
    % blah= 1;
else
    b11w_Lin= 0;
    b15w_Lin= 0;
    b51w_Lin= 0;
    b55w_Lin= 0;
end
B(1,1)= B(1,1)+b11w_Lin;
B(1,5)= B(1,5)+b15w_Lin;
B(5,1)= B(5,1)+b51w_Lin;
B(5,5)= B(5,5)+b55w_Lin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if strcmp(domainSwitch, 'time')
    %inputs for model:
    %linear: w tables, A table, B table, F table, C, I, Blin, numDOF
    %nonlinear: viscous damping, generator torque, wind forcing,
    %control pitch, Nonlin mooring lines
    
    %REQUIRED SIMULATE OUTPUTS:  [D_state, n_state, Stress_Tower, sigmaMean, X_RMS, V_RMS, array_Power, array_Power2, sigma_Nacelle, windPower]
    
    %output: time series of each DOF
    [D_state, n_state, Stress_Tower, sigmaMean, X_RMS, V_RMS, array_Power, array_Power2, sigma_Nacelle, windPower]=...
        simulateFWT_timeDomain_20170512(M, B, C, w_p, Hs,seaType,...
        wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable, ...
        A11Table, A15Table, A22Table, A24Table, A33Table, A24Table, A44Table, A15Table, A55Table, A66Table, ...
        B11Table, B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table,...
        rho_a, S, Rrotor, L_SWL_to_Nacelle, xDC, BetaMean, U, OmegaMean, Ngear, Idrivetrain,...
        qW_Matrix, Omega_Matrix, Ct_Matrix, Cq_Matrix, Cp_Matrix,...
        zzP, D, Cd,H,...
        angle_Kp, KpTable, angle_Ki, KiTable, omegaT, Tgen, ...
        materialParam, p);
    
    
else
    %% %%%%%%%%%%%%%%%% ITERATE UNTIL CONVERGED SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Borig= B;
    Corig= C; %B and C are reset to orig for every frequency
    while converge== 0 %set converge= 1 if responses accounting for nonlinearities have converged
        %if totally linear, don't iterate
        if (strcmp(mSwitch, 'lin') && strcmp(cSwitch, 'lin') && Cd== 0 && kzn== 0)
            converge= 1;
        end
        C= Corig;
        B= Borig;
        
        %% %%%% EQUIVALENT VISCOUS DAMPING COEFFICIENT %%%%%%%%%%%%%
        c11e= rho*Cd*(2/pi)^.5*sum(D.*sigmaQdot').*dzzP;
        %c15e= rho*Cd*(2/pi)^.5*sum((zzP-zCM).*D.*sigmaQdot').*dzzP;
        c15e= rho*Cd*(2/pi)^.5*sum(zzP.*D.*sigmaQdot').*dzzP;
        c51e= -c15e;
        %c55e= rho*Cd*(2/pi)^.5*sum((zzP-zCM).^2.*D.*sigmaQdot').*dzzP;
        c55e= -rho*Cd*(2/pi)^.5*sum(zzP.^2.*D.*sigmaQdot').*dzzP;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(cSwitch, 'Nonlin') && OmegaMean==0
            windPower=0;
        end
        
        %% %%%% EQUIVALENT NONLINEAR GENERATOR TORQUE AND WIND FORCING %%%%%%%%%%%%%%%
        if strcmp(cSwitch, 'Nonlin') && OmegaMean>0
            sO= Sigma(1,1)^.5;
            Omega_vector2= linspace(OmegaMean-4*sO, OmegaMean+4*sO,49);
            dOm2= Omega_vector2(2)-Omega_vector2(1);
            
            sQ= Sigma(2,2)^.5;
            qW_vector2= linspace(U-4*sQ, U+4*sQ,50);
            dqW2= qW_vector2(2)-qW_vector2(1);
            
            %     intOmega_vector2= linspace(4*sOmInt, 4*sOmInt,50);
            
            %     qDotwMax= U+sQ;
            if 1%(qDotwMax<11) %if beta kept constant
                mu= [OmegaMean, U];
                
                %sB= Sigma(3,3)^.5;
                %Beta_vector2= linspace(BetaMean-3*sO, BetaMean+3*sO,48); %BETA IS IN DEG
                %dB2= Beta_vector2(2)-Beta_vector2(1);
                
                [Omega_Matrix2, qW_Matrix2] = ndgrid(Omega_vector2, qW_vector2);
                
                fqwOm = multivariateNormal(mu, Sigma(1:2,1:2), Omega_Matrix2, qW_Matrix2); %order matters!
                fOm2= normpdf(Omega_vector2, OmegaMean, sO);
                %fBeta2= normpdf(Beta_vector2, BetaMean, sB); %BETA IS IN DEG
                
                Omega_vector2(Omega_vector2<0)= 0; %Adjustment for interpolation
                Omega_vector2(Omega_vector2>1.511)= 1.511;
                Omega_Matrix2(Omega_Matrix2<0)= 0; %Adjustment for interpolation
                Omega_Matrix2(Omega_Matrix2>1.511)= 1.511;
                qW_Matrix2(qW_Matrix2<min(min(qW_Matrix)))= min(min(qW_Matrix)); %Adjustment for interpolation
                qW_Matrix2(qW_Matrix2>max(max(qW_Matrix)))= max(max(qW_Matrix));
                
                %[qW_Matrix, dqW, Omega_Matrix, dOm, Ct_Matrix, Cq_Matrix, dCt_dq_Matrix, dCq_dq_Matrix, dCt_dOm_Matrix, dCq_dOm_Matrix]= myThrust_versus_Urel_Om_Matrix(beta);
                Ct_Matrix2= interp2(qW_Matrix, Omega_Matrix, Ct_Matrix, qW_Matrix2,Omega_Matrix2, 'cubic'); %interp2(X,Y,V,Xq,Yq,'cubic'); %X, Y, Xq, Yq are query grid.
                Cq_Matrix2= interp2(qW_Matrix, Omega_Matrix, Cq_Matrix, qW_Matrix2,Omega_Matrix2,'cubic');
                Cp_Matrix2= interp2(qW_Matrix, Omega_Matrix, Cp_Matrix, qW_Matrix2,Omega_Matrix2,'cubic');
                %       dCt_dq_Matrix2= interp2(qW_Matrix, Omega_Matrix, dCt_dq_Matrix, qW_Matrix2,Omega_Matrix2,'cubic');
                %       dCq_dq_Matrix2= interp2(qW_Matrix, Omega_Matrix, dCq_dq_Matrix, qW_Matrix2,Omega_Matrix2,'cubic');
                dCt_dOm_Matrix2= interp2(qW_Matrix, Omega_Matrix, dCt_dOm_Matrix, qW_Matrix2,Omega_Matrix2,'cubic');
                dCq_dOm_Matrix2= interp2(qW_Matrix, Omega_Matrix, dCq_dOm_Matrix, qW_Matrix2,Omega_Matrix2,'cubic');
                dTdO= interp1(omegaT, gradTgen, Omega_vector2);
                
                E_Ct_absQ= sum(sum( Ct_Matrix2.*abs(qW_Matrix2).*fqwOm  ))*dqW2*dOm2;
                %       E_absQ_Q_dCtdQ= sum(sum( dCt_dq_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOm  ))*dqW2*dOm2;
                E_Cq_absQ= sum(sum( Cq_Matrix2.*abs(qW_Matrix2).*fqwOm  ))*dqW2*dOm2;
                %       E_absQ_Q_dCqdQ= sum(sum( dCq_dq_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOm  ))*dqW2*dOm2;
                E_absQ_Q_dCtdO= sum(sum( dCt_dOm_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOm  ))*dqW2*dOm2;
                E_absQ_Q_dCqdO= sum(sum( dCq_dOm_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOm  ))*dqW2*dOm2;
                E_dTdO= sum(dTdO.*fOm2).*dOm2;
                %       E_absQ_Q_dCtdB=0; %Expected gradient of thrust coefficient wrt beta
                E_absQ_Q_dCqdB= 0; %Expected gradient of torque coefficient wrt beta
                E_dFpdOm= 0;
                E_dFidOm= 0;
                
                
                %Changes on April 24, 2017
                b1_1w= rho_a*S*E_Ct_absQ;
                b1_5w= L_SWL_to_Nacelle*b1_1w;
                b5_1w= L_SWL_to_Nacelle*b1_1w;
                b5_5w= L_SWL_to_Nacelle^2*b1_1w;
                b11_1w= rho_a*S*E_Cq_absQ;
                b11_5w= L_SWL_to_Nacelle*b11_1w;
                
                
                k1_11w= -.5*rho_a*S*E_absQ_Q_dCtdO;
                k1_13w= -.5*rho_a*S*E_absQ_Q_dCtdO;
                k5_11w= -L_SWL_to_Nacelle.*.5*rho_a*S*E_absQ_Q_dCtdO;
                k5_13w= -L_SWL_to_Nacelle.*.5*rho_a*S*E_absQ_Q_dCtdO;
                k11_11w= -.5*rho_a*S*Rrotor*E_absQ_Q_dCqdO - Ngear*E_dTdO;
                k11_13w= -.5*rho_a*S*Rrotor*E_absQ_Q_dCqdB;
                k13_11w= -E_dFpdOm;
                k13_13w= 0;
                L13_11w= -E_dFidOm;
                L13_13w= -E_dFidOm;
                
                
                %         k1_12w= -.5*rho_a*S*(E_Ct_absQ + E_absQ_Q_dCtdQ);
                %         k5_12w= L_SWL_to_Nacelle*k1_12w;
                %         k11_12= -.5*rho_a*S*Rrotor*(E_Cq_absQ + E_absQ_Q_dCqdQ);
                %         k1_11= -.5*rho_a*S*E_absQ_Q_dCtdO;
                %         k11_11= -.5*rho_a*S*E_absQ_Q_dCqdO + Ngear*E_dTdO;
                %
                %         %related to beta
                %         k1_13= -.5*rho_a*S*E_absQ_Q_dCtdB;
                %         k11_13= -.5*rho_a*S*E_absQ_Q_dCqdB;
                %         k5_11= L_SWL_to_Nacelle*k1_11;
                %         k5_13= L_SWL_to_Nacelle*k1_13;
                %         k13_11= 0; %-Ngear*E_Kp;
                %         L13_11= 0; %-Ngear*E_Ki;
                
                %Average Wind Power
                EqW3Cp= sum(sum( Cp_Matrix2.*qW_Matrix2.^3.*fqwOm  ))*dqW2*dOm2;
                windPower= .5*rho_a*S*EqW3Cp./10^6; %[MW]. Power= .5rhoA*S*E{Cp*qwDot^3}
                
            else %beta changes
                %       mu= [OmegaMean, U, BetaMean]; %x order is [Omega, qDot, beta]
                
                if Sigma(3,3)==0 %#ok<UNRCH>
                    Sigma(3,3)= .25;
                end
                sB= Sigma(3,3)^.5;
                Beta_vector2= linspace(BetaMean-3*sB, BetaMean+3*sB,48); %BETA IS IN DEG
                dB2= Beta_vector2(2)-Beta_vector2(1);
                
                [Omega_Matrix2, qW_Matrix2, Beta_Matrix2] = ndgrid(Omega_vector2, qW_vector2, Beta_vector2);
                %Each    row (x-value) is different qW.
                %Each column (y-value) is different Omega.
                %Each   page (z-value) is different beta.
                
                fqwOmB = multivariateNormal3(mu, Sigma, Omega_Matrix2, qW_Matrix2, Beta_Matrix2); %order matters!
                fOm2= normpdf(Omega_vector2, OmegaMean, sO);
                fBeta2= normpdf(Beta_vector2, BetaMean, sB); %BETA IS IN DEG
                
                
                Omega_vector2(Omega_vector2<0)= 0; %Adjustment for interpolation
                Omega_vector2(Omega_vector2>1.511)= 1.511;
                Omega_Matrix2(Omega_Matrix2<0)= 0; %Adjustment for interpolation
                Omega_Matrix2(Omega_Matrix2>1.511)= 1.511;
                qW_Matrix2(qW_Matrix2<min(min(qW_Matrix)))= min(min(qW_Matrix)); %Adjustment for interpolation
                qW_Matrix2(qW_Matrix2>max(max(qW_Matrix)))= max(max(qW_Matrix));
                
                %NOTE: Torque-rate limit of 15,000 Nm/s
                %Tgen2= interp1(omegaT, Tgen, Omega_vector2);
                
                %[qW_Matrix3, dqW3, Omega_Matrix3, dOm3, Ct_Matrix3, Cq_Matrix3, dCt_dq_Matrix3, dCq_dq_Matrix3, dCt_dOm_Matrix3, dCq_dOm_Matrix3, dCt_dB_Matrix3, dCq_dB_Matrix3]= myThrust_versus_Urel_Om_Matrix3;
                %NEED MESHGRID FORMAT
                Ct_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, Ct_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic'); %interp2(X,Y,V,Xq,Yq,'cubic'); %X, Y, Xq, Yq are query grid.
                Cq_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, Cq_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                Cp_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, Cp_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                dCt_dq_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, dCt_dq_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                dCq_dq_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, dCq_dq_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                dCt_dOm_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, dCt_dOm_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2,'cubic');
                dCq_dOm_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, dCq_dOm_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                dCt_db_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, dCt_dB_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                dCq_db_Matrix2= interp3(qW_Matrix3, Omega_Matrix3, Beta_Matrix3, dCq_dB_Matrix3, qW_Matrix2, Omega_Matrix2, Beta_Matrix2, 'cubic');
                
                dTdO= interp1(omegaT, gradTgen, Omega_vector2);
                %[angle_Kp, KpTable, angle_Ki, KiTable]= Jonkman_bladeGain;
                Kp2= interp1(angle_Kp, KpTable, Beta_vector2);
                Ki2= interp1(angle_Ki, KiTable, Beta_vector2);
                Ki2(isnan(Ki2))= KiTable(1);
                Kp2(isnan(Kp2))= KpTable(1);
                
                Ct_Matrix2(isnan(Ct_Matrix2))= 0;
                Cq_Matrix2(isnan(Cq_Matrix2))= 0;
                Cp_Matrix2(isnan(Cp_Matrix2))= 0;
                dCq_dOm_Matrix2(isnan(dCq_dOm_Matrix2))= 0;
                dCq_db_Matrix2(isnan(dCq_db_Matrix2))= 0;
                dCq_dq_Matrix2(isnan(dCq_dq_Matrix2))= 0;
                dCt_dOm_Matrix2(isnan(dCt_dOm_Matrix2))= 0;
                dCt_dq_Matrix2(isnan(dCt_dq_Matrix2))= 0;
                dCt_db_Matrix2(isnan(dCt_db_Matrix2))= 0;
                
                %SEPT 13
                %expectations are functions of qW, Omega, Beta
                E_Ct_absQ= sum(sum(sum( Ct_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2;
                E_absQ_Q_dCtdQ= sum(sum(sum( dCt_dq_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2;
                E_Cq_absQ= sum(sum(sum( Cq_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2;
                E_absQ_Q_dCqdQ= sum(sum(sum( dCq_dq_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2;
                E_absQ_Q_dCtdO= sum(sum(sum( dCt_dOm_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2;
                E_absQ_Q_dCqdO= sum(sum(sum( dCq_dOm_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2;
                E_dTdO= sum(dTdO.*fOm2).*dOm2;
                
                E_absQ_Q_dCtdB= sum(sum(sum( dCt_db_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2; %Expected gradient of thrust coefficient wrt beta
                E_absQ_Q_dCqdB= sum(sum(sum( dCq_db_Matrix2.*qW_Matrix2.*abs(qW_Matrix2).*fqwOmB  )))*dqW2*dOm2*dB2; %Expected gradient of torque coefficient wrt beta
                %         E_Kp= sum(Kp2.*fBeta2).*dB2;
                %         E_Ki= sum(Ki2.*fBeta2).*dB2;
                
                Kp= .00628./(1+(Beta_Matrix2.*pi./180./.11)); %Beta is in Deg
                Kp(Omega_Matrix2<122)= 0;
                Ki= .0008965./(1+(Beta_Matrix2.*pi./180./.11)); %Beta is in Deg
                Ki(Omega_Matrix2<122)= 0;
                dKp_dBeta= -0.0571./(9.09.*Beta_Matrix2+1).^2;
                dKi_dBeta= -0.00815/(9.09.*Beta_Matrix2+1).^2;
                
                %dFp_dOm= -Ngear.*Kp;
                OmTarget= 122.91;
                E_dFpdOm= -Ngear.*sum(sum(sum( Kp.* fqwOmB )   )    )*dqW2*dOm2*dB2; %Force of proportional feedback
                E_dFpdBeta= -Ngear.*sum(sum(sum( -dKp_dBeta.*(Omega_Matrix2-OmTarget).*fqwOmB  ) ) )*dqW2*dOm2*dB2; %Force of proportional feedback
                
                E_dFidOm= -Ngear.*sum(sum(sum( -Ki.* fqwOmB )   )    )*dqW2*dOm2*dB2;
                E_dFidBeta= -Ngear.*sum(sum(sum( -dKi_dBeta.*(Omega_Matrix2-OmTarget).*fqwOmB  ) ) )*dqW2*dOm2*dB2;
                
                %RAO coefficients related to dqW and Om DOFs
                %         k1_12w= -.5*rho_a*S*(E_Ct_absQ + E_absQ_Q_dCtdQ);
                %         k5_12w= L_SWL_to_Nacelle*k1_12w;
                %         k11_12= -.5*rho_a*S*Rrotor*(E_Cq_absQ + E_absQ_Q_dCqdQ);
                %         k1_11= -.5*rho_a*S*E_absQ_Q_dCtdO;
                %         k11_11= -.5*rho_a*S*E_absQ_Q_dCqdO + Ngear*E_dTdO;
                %
                %         %RAO coefficients related to beta DOFs
                %         k1_13= -.5*rho_a*S*E_absQ_Q_dCtdB;
                %         k11_13= -.5*rho_a*S*E_absQ_Q_dCqdB;
                %         k5_11= L_SWL_to_Nacelle*k1_11;
                %         k5_13= L_SWL_to_Nacelle*k1_13;
                %         k13_11= -Ngear*E_Kp;
                %         L13_11= -Ngear*E_Ki;
                
                
                %Changes on April 24, 2017
                b1_1w= rho_a*S*E_Ct_absQ;
                b1_5w= L_SWL_to_Nacelle*b1_1w;
                b5_1w= L_SWL_to_Nacelle*b1_1w;
                b5_5w= L_SWL_to_Nacelle^2*b1_1w;
                b11_1w= rho_a*S*E_Cq_absQ;
                b11_5w= L_SWL_to_Nacelle*b11_1w;
                
                
                k1_11w= -.5*rho_a*S*E_absQ_Q_dCtdO;
                k1_13w= -.5*rho_a*S*E_absQ_Q_dCtdO;
                k5_11w= -L_SWL_to_Nacelle.*.5*rho_a*S*E_absQ_Q_dCtdO;
                k5_13w= -L_SWL_to_Nacelle.*.5*rho_a*S*E_absQ_Q_dCtdO;
                k11_11w= -.5*rho_a*S*Rrotor*E_absQ_Q_dCqdO - Ngear*E_dTdO;
                k11_13w= -.5*rho_a*S*Rrotor*E_absQ_Q_dCqdB;
                k13_11w= -E_dFpdOm;
                k13_13w= -E_dFpdBeta;
                L13_11w= -E_dFidOm;
                L13_13w= -E_dFidBeta;
                
                %Average Wind Power
                EqW3Cp= sum(sum(sum( Cp_Matrix2.*qW_Matrix2.^3.*fqwOmB  )))*dqW2*dOm2*dB2;
                windPower= .5*rho_a*S*EqW3Cp./10^6; %[MW]. Power= .5rhoA*S*E{Cp*qwDot^3}
            end
            
            
        else %otherwise, Do not have control AND wind > 0. Ignore wind damping effects
            %         k1_12w= 0;
            %         k5_12w= 0;
            %         k11_12= 0;
            %         k1_11= 0;
            %         k1_13= 0;
            %         k11_13= 0;
            %         k5_11= 0;
            %         k5_13= 0;
            %         k11_11= 0;
            %         k13_11= 0;
            %         L13_11= 0;
            %         %windPower= 0;
            
            
            b1_1w= 0;
            b1_5w= 0;
            b5_1w= 0;
            b5_5w= 0;
            b11_1w= 0;
            b11_5w= 0;
            
            
            k1_11w= 0;
            k1_13w= 0;
            k5_11w= 0;
            k5_13w= 0;
            k11_11w= 0;
            k11_13w= 0;
            k13_11w= 0;
            k13_13w= 0;
            
            L13_11w= 0;
            L13_13w= 0;
            
        end
        
        % C(1,12)= C(1,12)+ k1_12w;
        % C(5,12)= C(5,12)+ k5_12w;
        % C(11,12)= C(11,12)+ k11_12;
        % C(1,11)= C(1,11) + k1_11;
        % C(1,13)= C(1,13) + k1_13;
        % C(5,11)= C(5,11) + k5_11;
        % C(5,13)= C(5,13) + k5_13;
        % C(11,11)= C(11,11) + k11_11;
        % C(11,13)= C(11,13) + k11_13;
        % C(13,11)= C(13,11) + k13_11;
        
        B(1,1)= B(1,1) + b1_1w;
        B(1,5)= B(1,5) + b1_5w;
        B(5,1)= B(5,1) + b5_1w;
        B(5,5)= B(5,5) + b5_5w;
        B(11,1)= B(11,1) + b11_1w;
        B(11,5)= B(11,5) + b11_5w;
        
        C(1,11)= C(1,11) + k1_11w;
        C(1,13)= C(1,13) + k1_13w;
        C(5,11)= C(5,11) + k5_11w;
        C(5,13)= C(5,13) + k5_13w;
        C(11,11)= C(11,11) + k11_11w;
        C(11,13)= C(11,13) + k11_13w;
        C(13,11)= C(13,11) + k13_11w;
        C(13,13)= C(13,13) + k13_13w;
        
        LL(13,11)= L13_11w;
        LL(13,13)= L13_13w;
        
        %Linearized wind damping moved to line above call for time domain simulation
        
        %% %%%%%%%%%%%%%%%%%%%%% NONLINEAR MOORING LINE %%%%%%%%%%%%%%%
        if strcmp(mSwitch, 'Nonlin')
            %Q is horizontal displacement at fair leads
            qq_Vector= linspace(qDC-4*sigmaQ, qDC+4*sigmaQ,200);
            dqq= qq_Vector(2)-qq_Vector(1);
            
            %sigmaQ, meanQ, q_table, f_q
            f_q= 1./((2.*pi).^.5.*sigmaQ).*exp( -(qq_Vector-qDC).^2./(2.*sigmaQ.^2) ); %Gaussian probability distribution
            dFdqIt= interp1(q_table, dFdq, qq_Vector);
            
            % dFdq= -4.118e4;
            E_dFdq= sum(dFdqIt.*f_q).*dqq;
            % E_dFdq= -4.118e4;
            E_dFdq(isnan(E_dFdq))= 0; %during first iteration, equals 0
            
            c11m= -E_dFdq;
            c15m= -Zfl.*E_dFdq;
            c55m= -Zfl.^2.*E_dFdq;
        else
            c11m= 0;
            c15m= 0;
            c55m= 0;
        end
        
        %% %%%% NONLINEAR TMDe effects
        if TMDe %djznEq
            if kjxn>0
                kjxnEq= 3*kjxn*Sigma4JTMDe.^2;
            else
                KjxnEq= 0;
            end
            if kjyn>0
                kjynEq= 3*kjyn*Sigma5JTMDe.^2;
            else
                KjynEq= 0;
            end
            for ii= 1:numTMDe
                if kzn>0
                    %          if strcmp(name,'mono')
                    %              K33nEq=
                    %          else
                    K33nEq(ii)= 3*kzn*SigmaWECheave(ii).^2; %#ok<AGROW>
                    %          end
                else
                    K33nEq(ii)= 0; %#ok<AGROW>
                end
            end
            
            %     if dj11nTMDe>0
            %         dj11nEq= dj11nTMDe*(2/pi)^.5/Sigma4JTMDe;
            %     else
            %         dj11nEq= 0;
            %     end
            %     if dj33nTMDe>0
            %         dj33nEq= dj33nTMDe*(2/pi)^.5/Sigma4JTMDe;
            %     else
            %         dj33nEq= 0;
            %     end
            
            if djxn>0
                djxnEq= djxn*(2/pi)^.5/Sigma4JTMDe;
            else
                djxnEq= 0;
            end
            if djyn>0
                djynEq= djxn*(2/pi)^.5/Sigma5JTMDe;
            else
                djynEq= 0;
            end
            if djzn>0
                djznEq= djzn*(2/pi)^.5/Sigma6JTMDe;
            else
                djznEq= 0;
            end
            
        end
        AorigII= A;
        BorigII= B;
        CorigII= C;
        %% ITERATE THROUGH FREQUENCIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for q= 1:length(wVector) %Need to iterate through w because each w has 2D matrices
            w= wVector(q);
            k= kVector(q);
            Vg= VgVector(q);
            B= BorigII;
            C= CorigII;
            A= AorigII;
            X= zeros(numDOF,1); %forcing vector
            X1TMDes= 0;
            X5TMDes= 0;
            % A= zeros(numDOF, numDOF)
            %     Cyl= 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %%% RAO MATRICES AT EACH FREQUENCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LL(13,11)= L13_11w;
            
            if useWpAB
                wUse= w_p;
            else
                wUse= w;
            end
            
            A(1,1)= A(1,1) + interp1(wTable, A11Table, wUse)    + A11Cyl;
            A(1,5)= A(1,5) + interp1(wTable, A15Table, wUse)    + A15Cyl;
            A(5,1)= A(1,5);
            A(2,2)= A(2,2) + interp1(wTable, A22Table, wUse)    +A22Cyl;
            A(2,4)= A(2,4) + interp1(wTable, A24Table, wUse)    +A24Cyl;
            A(4,2)= A(2,4);
            A(3,3)= A(3,3) + interp1(wTable, A33Table, wUse)    +A33Cyl;
            A(4,4)= A(4,4) + interp1(wTable, A44Table, wUse)    +A44Cyl;
            A(5,5)= A(5,5) + interp1(wTable, A55Table, wUse)    +A55Cyl;
            A(6,6)= A(6,6) + interp1(wTable, A66Table, wUse)    +A66Cyl;
            I= M+A;
            
            % B11Table= B11Table+ 1e5; %[Ns/m]
            % B22Table= B22Table+ 1e5;
            % B33Table= B33Table+ 1.3e5;
            % B66Table= B66Table+ 13e6;
            
            B(1,1)= B(1,1)+ interp1(wTable, B11Table, wUse)+  c11e;
            B(1,5)= B(1,5)+ interp1(wTable, B15Table, wUse)+  c15e;
            B(5,1)= B(5,1)+ interp1(wTable, B15Table, wUse)+  c51e;
            B(2,2)= B(2,2)+ interp1(wTable, B22Table, wUse);
            B(2,4)= B(2,4)+ interp1(wTable, B24Table, wUse);
            B(4,2)= B(4,2)+ interp1(wTable, B24Table, wUse);
            B(3,3)= B(3,3)+ interp1(wTable, B33Table, wUse);
            B(4,4)= B(4,4)+ interp1(wTable, B44Table, wUse);
            B(5,5)= B(5,5)+ interp1(wTable, B55Table, wUse)+  c55e;
            B(6,6)= B(6,6)+ interp1(wTable, B66Table, wUse);
            
            %Additional damping from sphere is defined below in the forcing section
            
            if TMDe
                for ii=1:numTMDe
                    %     xwInd= 14+5*ii;
                    %     ywInd= 15+5*ii;
                    zwInd= 11+3*ii;
                    %     zColInd= 17+5*ii;
                    %     PInd= 18+5*ii;
                    %
                    %     C(1,1)=      C(1,1) + K11nEq;
                    %     C(1,5)=      C(1,5)+ K11nEq*zTMDeAttatch;
                    %     C(1,xwInd)= C(1,xwInd)- K11nEq;
                    %     C(5,1)=      C(5,1) + K11nEq*zTMDeAttatch;
                    %     C(5,5)=      C(5,5) + K11nEq*zTMDeAttatch^2;
                    %     C(5,xwInd)= C(5,14+TMD)+ K11nEq*zTMDeAttatch;
                    %     C(xwInd, 1)= C(14+TMD, 1) -K11nEq;
                    %     C(xwInd, 5)= C(14+TMD, 5) -K11nEq*zTMDeAttatch;
                    %     C(xwInd, xwInd)= C(14+TMD, 14+TMD) + K11nEq;
                    %
                    %     C(3,3)=      C(3,3)      + K33nEq;
                    %     C(3,15+TMD)= C(3,15+TMD) -K33nEq;
                    %     C(15+TMD, 3)= C(15+TMD, 3) -K33nEq;
                    C(zwInd, zwInd)= C(zwInd, zwInd) + K33nEq(ii);
                    %
                    %     B(1,1)=      B(1,1) +dj11nEq;
                    %     B(1,5)=      B(1,5)+dj11nEq*zTMDeAttatch;
                    %     B(1,14+TMD)= B(1,14+TMD)-dj11nEq;
                    %
                    %     B(5,1)=      B(5,1) +dj11nEq*zTMDeAttatch;
                    %     B(5,5)=      B(5,5) +dj11nEq*zTMDeAttatch^2;
                    %     B(5,14+TMD)= B(5,14+TMD)+dj11nEq*zTMDeAttatch;
                    %
                    %     B(14+TMD, 1)= B(14+TMD, 1) -dj11nEq;
                    %     B(14+TMD, 5)= B(14+TMD, 5) -dj11nEq*zTMDeAttatch;
                    %     B(14+TMD, 14+TMD)= B(14+TMD, 14+TMD) +dj11nEq;
                    %
                    %
                    %     B(3,3)=      B(3,3)      + dj33nEq;
                    %     B(3,15+TMD)= B(3,15+TMD) -dj33nEq;
                    %
                    %     B(15+TMD, 3)= B(15+TMD, 3) -dj33nEq;
                    %     B(15+TMD, 15+TMD)= B(15+TMD, 15+TMD) + dj33nEq;
                end
            end
            
            
            %SPRING
            C(1,1)= C(1,1) + c11m;
            C(1,5)= C(1,5) + c15m;
            C(5,1)= C(5,1) + c15m;
            C(5,5)= C(5,5) + c55m;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %%% FORCING AT EACH FREQUENCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalized by wave amplitude
            
            %%% Equivalent Viscous forcing on platform
            % TF_Uwater= -g.*k./w.*cosh(k.*(zzP+H))./cosh(k.*H);
            X1e= 0;%rho*Cd*(2/pi)^.5*sum(D.*sigmaQdot'.*TF_Uwater).*dzzP;
            X5e= 0;%-rho*Cd*(2/pi)^.5*sum(zzP.*D.*sigmaQdot'.*TF_Uwater).*dzzP;
            
            
            % %%% Wells in spar forcing
            % if WellsinSpar ==1
            %     TF_Vwater= g.*k./w.*cosh(k.*(z_WellsinSpar+H))./cosh(k.*H);
            %     TF_Vwater2= g.*k./w.*cosh(k.*(z_WellsinSpar2+H))./cosh(k.*H);
            %
            %     Xwells1= +d_WellsinSpar.*TF_Vwater        + d_WellsinSpar2.*TF_Vwater2;
            %     %Xwells5= X(5) +d_WellsinSpar.*TF_Vwater.*L_pw  + d_WellsinSpar2.*TF_Vwater2.*L_pw2; %Lpw, z_A11 are negative values
            %     Xwells5= +d_WellsinSpar.*TF_Vwater.*z_WellsinSpar  + d_WellsinSpar2.*TF_Vwater2.*z_WellsinSpar2;
            % else
            %     Xwells1= 0;
            %     Xwells5= 0;
            % end
            
            %%% Cylinder forcing- Froude Krylov
            if Cyl
                TF_A1water= 1i.*g.*k.*cosh(k.*(zCyl+H))./cosh(k.*H);
                % TF_A1waterTMDe= 1i.*g.*k.*cosh(k.*(zS+H))./cosh(k.*H).*exp(1i.*k*xS);
                TF_A3water= -g.*k.*sinh(k.*(zCyl+H))./cosh(k.*H);
                %TF_A3waterTMDe= -g.*k.*sinh(k.*(zH+H))./cosh(k.*H).*exp(1i.*k*xH);
                X1Cyl= -(rho*VCylSub+A11Cyl).*TF_A1water;
                X3Cyl= -(rho*VCylSub+A33Cyl).*TF_A3water;
                X5Cyl= zCyl*X1Cyl; %zSphere is negative
                
                %Cylinder hydrodynamic damping
                B11Cyl= k/(8*rho*g*Vg)*abs(X1Cyl)^2;
                B33Cyl= k/(4*rho*g*Vg)*abs(X3Cyl)^2;
            elseif usePlate
                % Override cylinder parameters with heave plate parameters
                
                % Assume small volume
                
                TF_A1water= 1i.*g.*k.*cosh(k.*(heavePlate.Lfz+H))./cosh(k.*H);
                % TF_A1waterTMDe= 1i.*g.*k.*cosh(k.*(zS+H))./cosh(k.*H).*exp(1i.*k*xS);
                TF_A3water= -g.*k.*sinh(k.*(heavePlate.Lfz+H))./cosh(k.*H);
                %TF_A3waterTMDe= -g.*k.*sinh(k.*(zH+H))./cosh(k.*H).*exp(1i.*k*xH);
                X1Cyl= 0; %-(rho*VCylSub+A11Cyl).*TF_A1water;
                X3Cyl= -(heavePlate.A33TMDe).*TF_A3water;
                X5Cyl= 0;%heavePlate.Lfz*X1Cyl; %zSphere is negative
                
                %Cylinder hydrodynamic damping
                B11Cyl= k/(8*rho*g*Vg)*abs(X1Cyl)^2;
                B33Cyl= k/(4*rho*g*Vg)*abs(X3Cyl)^2;
            else
                X1Cyl= 0;
                X3Cyl= 0;
                X5Cyl= 0;
                B11Cyl= 0;
                B33Cyl= 0;
            end
            
            B(1,1)= B(1,1) + B11Cyl;
            B(2,2)= B(2,2) + B11Cyl;
            B(3,3)= B(3,3) + B33Cyl;
            B(1,5)= B(1,5) + zCyl*B11Cyl;
            B(5,1)= B(5,1) + zCyl*B11Cyl;
            B(2,4)= B(2,4) + -zCyl*B11Cyl;
            B(4,2)= B(4,2) + -zCyl*B11Cyl;
            B(5,5)= B(5,5) + zCyl^2*B11Cyl;
            B(4,4)= B(4,4) + zCyl^2*B11Cyl;
            
            
            % external tuned mass damper
            if TMDe
                for ii=1:numTMDe
                    
                    zS= zSVector; % identical for all array components
                    if length(zS)>1
                        dzS= zS(2) - zS(1);
                    else
                        dzS= 1;
                    end
                    xS= xSVector(:,ii); % different for each component
                    zH= zHVector;
                    xH= xHVector(:,ii); % diff
                    dVS= dVSVector;
                    dAS= dASVector;
                    dVH= dVHVector;
                    dAH= dAHVector;
                    
                    zWellsForceS= zWellsForceSVector(:,ii);
                    xWellsForceS= xWellsForceSVector(:,ii);
                    xWellsForceH= xWellsForceHVector(:,ii);
                    zWellsForceH= zWellsForceHVector(:,ii);
                    
                    %%% External TMD: Linear forcing: incident waves and Wells turbine),
                    TF_A1waterTMDe= 1i.*g.*k.*cosh(k.*(zS+H))./cosh(k.*H).*exp(1i.*k*xS);
                    TF_A3waterTMDe= -g.*k.*sinh(k.*(zH+H))./cosh(k.*H).*exp(1i.*k*xH);
                    TF_A3waterCol= -g.*k.*sinh(k.*(zWellsForceH+H))./cosh(k.*H).*exp(1i.*k*xWellsForceH);
                    dX1TMDe= -(rho.*dVS+dAS).*TF_A1waterTMDe; % Incremental Surge forcing
                    dX3TMDe= -(rho.*dVH+dAH).*TF_A3waterTMDe;
                    
                    X1TMDe= sum(dX1TMDe); %-(rho*VsubTMDe+A11TMDe).*TF_A1waterTMDe;
                    %X3TMDe= sum(dX3TMDe); %-(rho*VsubTMDe+A33TMDe).*TF_A3waterTMDe;
                    X3TMDe= C33TMDe; % Forcing is rho*g*A*a. June 30, 2019
                    X3Col= -(ICol).*TF_A3waterCol;
                    X5TMDe= sum(dX1TMDe.*zS);
                    
                    % TMDe
                    % Strip theory for B11, B15, B55
                    % 3D heave for B33's.
                    % Account for coordinate transformation here
                    B11TMDe= sum( abs( dX1TMDe ).^2 ) * dzS /(2*rho*g*Vg);
                    % can ignore B15 transformation term due to array symmetry?
                    % Used to check hinge force
                    B15TMDe= sum( abs( dX1TMDe ).^2 .* zS ) * dzS /(2*rho*g*Vg);
                    B55TMDe= sum( abs( dX1TMDe ).^2 .* zS.^2 ) * dzS /(2*rho*g*Vg) + B11TMDe*x_junction_WEC(ii).^2;
                    B55TMDe= sum( abs( dX1TMDe ).^2 .* zS.^2 ) * dzS /(2*rho*g*Vg) + B11TMDe*x_junction_WEC(ii).^2;
                    % B55~ z^2*B11- torque due to translation due to rotation
                    
                    B33TMDe= k/(4*rho*g*Vg)*abs(X3TMDe)^2; % Why differs by factor of 2 from surge?
                    B33Col= 25*k/(4*rho*g*Vg)*abs(X3Col)^2/ACol;
                    B33Col(isnan(B33Col))= 0; %Column not used in simulation
                    
                    %%% Final damping elements:
                    
                    % add rigid properties to FWT platform:
                    % Lines above account for coordinate transformation
                    B(1:6,1:6)= B(1:6,1:6) + ...
                        [B11TMDe     0    0       0   B15TMDe   0
                        0    B11TMDe   0    B15TMDe    0     0
                        0        0     0       0       0     0
                        0    B15TMDe   0    B55TMDe    0     0
                        B15TMDe     0     0       0    B55TMDe  0
                        0        0     0       0       0     0];
                    
                    % WEC's free properties:
                    zwInd= 11+3*ii;
                    zColInd= 12+3*ii;
                    PInd= 13+3*ii;
                    
                    B(zwInd, zwInd)= B(zwInd, zwInd)+ B33TMDe;
                    B(zColInd, zColInd)= B(zColInd, zColInd)+ B33Col;
                    B(zColInd, zwInd)= B(zColInd, zwInd)+ B33Col;
                    
                    %%% Exposed damper forcing
                    TF_Vxwater= g.*k.*w.*cosh(k.*(zWellsForceS+H))./cosh(k.*H).*exp(1i.*k*min(xWellsForceS));
                    X1absTMDe= dWellsx.*TF_Vxwater;
                    
                    TF_Vzwater= 1i.*g.*k./w.*sinh(k.*(zWellsForceH+H))./cosh(k.*H).*exp(1i.*k*min(xWellsForceH));
                    X3absTMDe= dWellsz.*TF_Vzwater;
                    
                    %%% Final force vector
                    X1TMDes= X1TMDes + X1TMDe + X1absTMDe;
                    X5TMDes= X5TMDes + X5TMDe;
                    % ignore array roll excitation due to symmetry
                    
                    X(zwInd)= X3TMDe + X3absTMDe;
                    if ACol>0
                        X(zColInd)= X3Col/ACol;
                    else
                        X(zColInd)= 0;
                    end
                    
                    
                    %VISCOUS FORCING
                    %TF_Uwater= g.*k./w.*cosh(k.*(zzP+H))./cosh(k.*H);
                    % X1e= rho*Cd*(2/pi)^.5*sum(D.*sigmaQdot'.*TF_Uwater).*dzzP;
                    % X5e= -rho*Cd*(2/pi)^.5*sum(zzP.*D.*sigmaQdot'.*TF_Uwater).*dzzP;
                    
                    
                end
                
            end
            
            %Wells turbine in spar
            TF_VxWells= g.*k.*w.*cosh(k.*(zWells+H))./cosh(k.*H);
            X1Wells= dWells.*TF_VxWells;
            
            
            %Final platform excitation vector
            X(1)= X(1) + F1(q)  + X1Cyl  + X1e   + X1Wells         + X1TMDes;
            X(3)= X(3) + F3(q)  + X3Cyl;
            X(4)= X(4);
            X(5)= X(5) + F5(q)  + X5Cyl + X5e    + X1Wells.*zWells + X5TMDes;
            X(6)= X(6);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% SOLVE LINEAR RAO AT EACH FREQUENCY %%%%%%%%%%%%%%%%%%%%%%%
            RAO_denom= -w^2.*I + 1i.*w.*B + C + LL./(1i.*w); %denominator of RAO (numerator is force per amplitude)
            
            %Wiener-Khinchine: Variance= integral:(S+)|H|^2 dw:
            alpha= pinv(RAO_denom);
            
            TF= RAO_denom\X; %alpha*X;
            S_response(:,q)= (abs( TF )).^2.*Su(q); %Power spectral density of response, for this sea state.
            TF_terms(:,q)= TF; %response per wave amplitude    runTypee
            X_terms(:,q)= X;
            
            %% Check that link heave force is negligible
            if strcmp(runType, 'compute_hinges')
                % 1. Calculate load in bar that is in Y-Z plane with FWT
                % WEC EOM's: In Coordinate system centered at FWT.
                % X1TMDe - X1Hinge = -w^2*(m+A11)X1 -w^2*A15*X5 + iw*B11*X1 + iw*B15*X5
                % X5TMDe - X5Hinge = -w^2*(I55+A55)X1 -w^2*A15*X1 + iw*B55*X5 + iw*B15*X1
                % FWT EOM's:
                
                F1bar_S(q)= abs((  -X1TMDe...
                    + ( -w^2*(mTMDe + A11TMDe(1)) + 1i*B11TMDe)*TF(1)...
                    + ( -w^2*(        A15TMDe(1)) + 1i*B15TMDe)*TF(5) )).^2.*Su(q);
                M5bar_S(q)= abs((  -X5TMDe...
                    + ( -w^2*(I55TMDe + A55TMDe(1)) + 1i*B15TMDe + C55TMDe)*TF(5)...
                    + ( -w^2*(          A51TMDe(1)) + 1i*B15TMDe)*TF(1) )).^2.*Su(q);
            end
            %% Matrix spectral density
            Sin= diag(X)*Su(q)*ctranspose(diag(X)); %matrix of forcing at this frequency. X is forcing normalized by amplitude
            S_responseMat= alpha*Sin*ctranspose(alpha); %MATCHES S_reponse!
            
            %Spectral outputs used for covariances for stat lin of turbine control
            SigmaTemp(:,q)= [S_responseMat(11,11); S_responseMat(11,12); S_responseMat(11,13); S_responseMat(12,12); S_responseMat(12,13); S_responseMat(13,13)];
            
            %Stress spectral density
            TF_Stress_surge(q)= (252*TF(7) -1.2081e4*TF(9)).*10^6; %Correction on Nov. 18, 2016. WHEN JONKMAN EIGENSHAPES [Pa]
            % TF_Stress_surge(q)= (252*TF(7)).*10^6; %MAY 11
            
            
        end
        %% SUM SOLUTIONS AT EACH FREQUENCY %%%%%%%%%%%%%%%%%%%%%%%
        SigmaTempSum= abs(sum(SigmaTemp, 2)*dw);
        Sigma= [SigmaTempSum(1) SigmaTempSum(2) SigmaTempSum(3);
            SigmaTempSum(2) SigmaTempSum(4) SigmaTempSum(5);
            SigmaTempSum(3) SigmaTempSum(5) SigmaTempSum(6)]; %covariance matrix for Omega, qDot, beta
        %Covariance matrix of
        
        %Require true for positive definite: all(eig(Sigma))>0. chol(Sigma)
        
        %Calculate sigmaQdot for viscous forcing and sigmaQ for mooring line forcing
        kMatrix= repmat(kVector,length(zzP),1); %kVector is row vector
        wMatrix= repmat(wVector,length(zzP),1);
        zzPMatrix= repmat(zzP',1,length(wVector));
        TFsurgeMatrix= repmat(TF_terms(1,:),length(zzP),1); %platform surge
        TFpitchMatrix= repmat(TF_terms(5,:),length(zzP),1); %platform pitch
        SuMatrix= repmat(Su, length(zzP), 1);
        func= cosh(kMatrix.*(zzPMatrix+H))./cosh(kMatrix.*H);
        func(isnan(func))= 1; %limit of cosh
        TF_Qdot= -kMatrix.*g./wMatrix.*func - 1i.*wMatrix.*TFsurgeMatrix- 1i.*zzPMatrix.*wMatrix.*TFpitchMatrix; %relative water velocity along platform
        TF_Q= TF_terms(1,:) + Zfl.*TF_terms(5,:); %horizontal displacement at fairleads
        sigmaQdot= (sum(abs(TF_Qdot).^2.*SuMatrix, 2).*dw).^.5; %column vector of st Dev in relative velocity along spar length [m/s]
        sigmaQ= (sum(abs(TF_Q).^2.*Su,2).*dw).^.5; %displacement at mooring line
        
        %TMDe-FWT junction spring
        %6 platform + 4 tower + Omega(11) + qDotw(12) + deltaBeta(13) + xTMD(14) + xTMDe(15) + xTMDe(16) + 2*number of OWC's (kinematic realationship)
        if TMDe
            for ii=1:numTMDe
                %         xwInd= 14+5*ii;
                %         ywInd= 15+5*ii;
                zwInd= 11+3*ii;
                %         zColInd= 17+5*ii;
                %         PInd= 18+5*ii;
                
                
                Lxl= x_FWT_to_WEC_vec(ii);
                Lyl= y_FWT_to_WEC_vec(ii);
                Lzl= z_FWT_to_WEC_vec; % same for all array components
                
                %     Dth= [dJx 0 0; 0 dJy 0; 0 0 dJz];
                JxtoTh= [0 -1/Lzl 1/Lyl; 1/Lzl 0 -1/Lxl; -1/Lyl 1/Lxl 0]; %theta= Jx * x
                %Lz= 0 and Ly= 0 mean no matter how move, there will never be an angle th_4
                JxtoTh(isinf(JxtoTh))= 0;
                Arel= [-1,0,0,0,0,0,1,0,0;0,-1,0,0,0,0,0,1,0;0,0,-1,0,0,0,0,0,1]; %Relative translations of link
                Arot= [0,0,0,1,0,0,0,0,0; 0,0,0,0,1,0,0,0,0; 0,0,0,0,0,1,0,0,0];
                
                thetaRel= (JxtoTh*Arel+Arot);
                
                %     B(4:6,1:6)= B(4:6,1:6)-TauRot(1:3,1:6);
                %     B(4:6,xwInd:zwInd)= B(4:6,xwInd:zwInd)-TauRot(1:3,7:9);
                %     B(xwInd:zwInd,1:6)= B(xwInd:zwInd,1:6)-FwecRot(1:3,1:6);
                %     B(xwInd:zwInd,xwInd:zwInd)= B(xwInd:zwInd,xwInd:zwInd)+FwecRot(1:3,7:9);
                
                %         TF4J= thetaRel(1,1).*TF_terms(1,:) + thetaRel(1,2).*TF_terms(1,:) + thetaRel(1,3).*TF_terms(3,:)+...
                %               thetaRel(1,4).*TF_terms(4,:) + thetaRel(1,5).*TF_terms(5,:) + thetaRel(1,6).*TF_terms(6,:)+...
                %               thetaRel(1,7).*TF_terms(xwInd,:) + thetaRel(1,8).*TF_terms(ywInd,:) + thetaRel(1,9).*TF_terms(zwInd,:);
                %
                %         TF5J= thetaRel(2,1).*TF_terms(1,:) + thetaRel(2,2).*TF_terms(2,:) + thetaRel(2,3).*TF_terms(3,:)+...
                %               thetaRel(2,4).*TF_terms(4,:) + thetaRel(2,5).*TF_terms(5,:) + thetaRel(2,6).*TF_terms(6,:)+...
                %               thetaRel(2,7).*TF_terms(xwInd,:) + thetaRel(2,8).*TF_terms(ywInd,:) + thetaRel(2,9).*TF_terms(zwInd,:);
                %
                %         TF6J= thetaRel(3,1).*TF_terms(3,:) + thetaRel(3,2).*TF_terms(3,:) + thetaRel(3,3).*TF_terms(3,:)+...
                %               thetaRel(3,4).*TF_terms(4,:) + thetaRel(3,5).*TF_terms(5,:) + thetaRel(3,6).*TF_terms(6,:)+...
                %               thetaRel(3,7).*TF_terms(xwInd,:) + thetaRel(3,8).*TF_terms(ywInd,:) + thetaRel(3,9).*TF_terms(zwInd,:);
                %         %TFZr= TF_terms(zwInd,:) - TF_terms(3,:);
                %
                %     Sigma4JTMDe= (sum(abs(TF4J).^2.*Su,2).*dw).^.5;
                %     Sigma5JTMDe= (sum(abs(TF5J).^2.*Su,2).*dw).^.5;%relative horizontal displacement between TMDe and platform
                %     Sigma6JTMDe= (sum(abs(TF6J).^2.*Su,2).*dw).^.5;    %SigmaZrTMDe= (sum(abs(TFZr).^2.*Su,2).*dw).^.5; %relative vertical displacement between TMDe and platform
                
                SigmaWECheave(ii)= (sum(abs(TF_terms(zwInd,:)).^2.*Su,2).*dw).^.5;
                
                %TFXpr= TFXr.*wVector;
                %         TFZpr= 0; %TFZr.*wVector;
                %SigmaXrPTMDe= (sum(abs(TFXpr).^2.*Su,2).*dw).^.5;%relative horizontal velocity between TMDe and platform
                %Sigma4JTMDe= (sum(abs(TF4J).^2.*Su,2).*dw).^.5; %relative vertical velocity between TMDe and platform
            end
        end
        
        %check convergence
        X_RMS_old_old= X_RMS_old;
        X_RMS_old= X_RMS;
        X_RMS= (sum(S_response, 2).*dw).^.5; %RMS of total motion. All frequency motions are independent because linear
        
        % Diffs1= abs((X_RMS(1)-X_RMS_old(1))/X_RMS_old(1));
        % Diffs5= abs((X_RMS(5)-X_RMS_old(5))/X_RMS_old(5));
        Diffs= abs((X_RMS-X_RMS_old)./X_RMS_old);
        diffs_oldOld= abs((X_RMS-X_RMS_old_old)./X_RMS_old_old); %means 2 solutions?
        maxDiffs= max(Diffs);
        maxDiffsOld= max(diffs_oldOld);
        %Check if responses have converged to within 5%
        if ( maxDiffs >.003 && it<80)
            %converge= 0;
        elseif( maxDiffs >.003 && it<80 )
            %converge= 0;
        else
            converge= 1;
        end
        it= it+1 %#ok<NOPRT>
        if it>30
            %         blah= 1;
            converge= 1;
        end
        if it>7 && maxDiffsOld<.003
            X_RMS= max(X_RMS, X_RMS_old);%Use the larger response
            converge= 1;
        end
        
    end
    
    
    %%
    %! S_response NOT already multiplied by dw
    
    %wMatrix_Col= repmat(wVector,N,1); %need to adjust when indep OWC. wMatrix_Col is used for OWC power
    wMatrix= repmat(wVector, numDOF, 1); %wMatrix is used for finding the RMS velocity
    
    Sv_response= wMatrix.^2.*S_response; %Power spectral density of velocities
    V_RMS= (sum(Sv_response, 2).*dw).^.5; %Total RMS velocity
    
    %  Sa_response= wMatrix.^4.*S_response; %A^2= w^4*X^2. Velocity of platform CM and nacelle bending
    % A_RMS= (sum(Sa_response, 2).*dw).^.5; %Total RMS acceleration
    
    Stress_Tower= (sum(abs(TF_Stress_surge).^2.*Su).*dw)^.5; %RMS stress. ONLY WORKS FOR NO SWAY. [Pa]
    Stress_Tower_Dot= (sum(abs(wVector.*TF_Stress_surge).^2.*Su).*dw)^.5;
    %Wiener-Khinchine: Variance= integral:(S+)|H|^2 dw:
    
    Fthrust= interp1(UbuoyTable, FthrustTable, U);
    sigmaMean= 90.*Fthrust; %march 24, 2016.   252/(1.6962e6).*Fthrust; %[Pa]
    
    %Stress in the tower in surge. For Constant Wind.
    if ~RAOrun %So don't need to declare tower material properties
        [D_state, n_state]= calculate_Fatigue(Stress_Tower, Stress_Tower_Dot, sigmaMean, p, seaType, materialParam); %damage to tower from this state. ASSUME SURGE BENDING ONLY FOR NOW
    else
        D_state= 0;
        n_state= 0;
    end
    
    %% Nacelle motion output
    TF_Nacelle= TF_terms(7,:) + TF_terms(9,:);
    S_Nacelle= (abs(TF_Nacelle)).^2.*Su; %S_response(:,q)= (abs( TF )).^2.*Su(q);
    sigma_Nacelle= (sum(S_Nacelle).*dw).^.5; %[m]
    
    
    
    %% Power output
    if TMDe
        Power_TMDe= 0;
        for ii= 1:numTMDe
            
            %         xwInd= 14+5*ii;
            %         ywInd= 15+5*ii;
            zwInd= 11+3*ii;
            %         zColInd= 17+5*ii;
            PInd= 13+3*ii;
            
            Lxl= x_FWT_to_WEC_vec(ii);
            Lyl= y_FWT_to_WEC_vec(ii);
            Lzl= z_FWT_to_WEC_vec;
            
            Lwz= z_junction_WEC;
            
            %         Dth= [dJx 0 0; 0 dJy 0; 0 0 dJz];
            JxtoTh= [0 -1/Lzl 1/Lyl; 1/Lzl 0 -1/Lxl; -1/Lyl 1/Lxl 0]; %theta= Jx * x
            %Lz= 0 and Ly= 0 mean no matter how move, there will never be an angle th_4
            JxtoTh(isinf(JxtoTh))= 0;
            Arel= [-1,0,0,0,0,0,1,0,0;0,-1,0,0,0,0,0,1,0;0,0,-1,0,0,0,0,0,1]; %Relative translations of link
            Arot= [0,0,0,1,0,0,0,0,0; 0,0,0,0,1,0,0,0,0; 0,0,0,0,0,1,0,0,0];
            
            thetaRel= (JxtoTh*Arel+Arot);
            
            S_Pc= abs(TF_terms(PInd,:)).^2.*Su;
            sigmaPC= (sum(S_Pc).*dw).^.5;
            
            %Translations
            TF_X11TMDe= TF_terms(1,:); % about SWL
            TF_V11TMDe= 1i.*wVector.*TF_X11TMDe;
            S_V11TMDe= abs(TF_V11TMDe).^2.*Su;
            %S_X11TMDe= abs(TF_X11TMDe).^2.*Su;
            Sigma_V11= (sum(S_V11TMDe).*dw).^.5;
            
            TF_X33TMDe= TF_terms(zwInd,:) -TF_terms(3,:);
            TF_V33TMDe= 1i.*wVector.*TF_X33TMDe;
            S_V33TMDe= abs(TF_V33TMDe).^2.*Su;
            %S_X33TMDe= abs(TF_X33TMDe).^2.*Su;
            Sigma_V33= (sum(S_V33TMDe).*dw).^.5;
            V33_WEC_wrt_Platform= Sigma_V33; % JK May 23, 2021
            
            TF_Vxwater= g.*kVector./wVector.*cosh(kVector.*(Lwz+H))./cosh(kVector.*H).*exp(1i.*kVector*min(xWellsForceS));
            %TF_Vabszwater= 1i.*g.*kVector./wVector.*sinh(kVector.*(zTMDeForce+H))./cosh(kVector.*H);
            
            % Power spectrum of relative horizontal velocity between wave and WEC
            %          S_absV11TMDe= abs(TF_Vxwater - 1i.*wVector.*TF_terms(xwInd,:)).^2.*Su;
            %          SigmaAbsV11= (sum(S_absV11TMDe).*dw).^.5;
            
            TF_Vzwater= 1i.*g.*kVector./wVector.*sinh(kVector.*(Lwz+H))./cosh(kVector.*H).*exp(1i.*kVector*min(xWellsForceH));
            S_absV33TMDe= abs(TF_Vzwater - 1i.*wVector.*TF_terms(zwInd,:)).^2.*Su;
            SigmaAbsV33= (sum(S_absV33TMDe).*dw).^.5;
            
            if kWells==0
                kWellsTemp= 1;
            else
                kWellsTemp= kWells;
            end
            %         Power_TMDe= Power_TMDe +  etaJ.*( (djxnEq+dJx).*sum(Sigma4JTMDe.^2)+ (djynEq+dJy).*sum(Sigma5JTMDe.^2)+...
            %                                 (djznEq+dJz).*sum(Sigma6JTMDe)) + etaW*(dWellsx.*SigmaAbsV11.^2 + dWellsz.*SigmaAbsV33.^2 + sigmaPC.^2./kWellsTemp) +...
            %                                 +etaJ.*(d1In.*Sigma_V11.^2 + d3In.*Sigma_V33.^2); %[W]
            
            % Power in W
            Power_TMDe= Power_TMDe +  ... %etaW*dWellsx.*SigmaAbsV11.^2 + ... Need to adjust coordinates to use this
                etaW*dWellsz.*SigmaAbsV33.^2 + ...
                etaW*sigmaPC.^2./kWellsTemp + ...
                etaJ.*(dJx.*Sigma_V11.^2 + dJz.*Sigma_V33.^2); %[W]
            
        end
    else
        Power_TMDe= 0;
        V33_WEC_wrt_Platform= 0;
        %     XrelTMDe= 0;
        %     ZrelTMDe= 0;
    end
    
    
    % Wells turbine
    TF_VxWells= g.*k.*wVector.*cosh(kVector.*(zWells+H))./cosh(kVector.*H);
    S_absWells= abs(TF_VxWells - 1i.*wVector.*TF_terms(1,:) - 1i.*wVector.*TF_terms(5,:).*zWells).^2.*Su;
    SigmaWells= (sum(S_absWells).*dw).^.5;
    etaW= .6;
    PWells= etaW*dWells*SigmaWells^2;
    
    
    %simulation output
    array_Power= array_Power   +  Power_TMDe  + PWells; %[W]
    array_Power2= Power_TMDe; %power of just damper 2
    
    
    % if OmegaMean==0
    %     windPower= 0;
    % elseif ~(strcmp(cSwitch, 'Nonlin') && OmegaMean>0)%~strcmp(cSwitch, 'Nonlin') %if linear
    %     sQ= Sigma(2,2);%standard deviation in qDot
    %    windPower= (.5*rho_a*S*Cp*U + 3*rho_a*S*Cp*U*sQ/2)./10^6;
    % end
    
    %% Hinge calculations
    if strcmp(runType, 'compute_hinges')
        M5bar_amp= (sum(M5bar_S).*dw).^.5 .* 2^.5;
        F1bar_amp= (sum(F1bar_S).*dw).^.5 .* 2^.5;
        
        hinge_stress_amp= zeros(100, 20);
        r_out_vec= linspace(0.2,2,100);
        r_in_frac_vec= linspace(0.9, 0.99, 20);
        for i= 1:length(r_out_vec)
            r_out= r_out_vec(i);
            for j= 1:length(r_in_frac_vec)
                r_in= r_in_frac_vec(j).*r_out;
                J_hinge= pi*(r_out^4-r_in^4)/2; % torsional rigidity
                I_hinge= J_hinge/2;
                array_radius= ( x_FWT_to_WEC_vec(1).^2 +y_FWT_to_WEC_vec(1).^2 ).^.5;
                
                % Von Mises Stress
                % Add shear stress and axial force stress?
                hinge_stress_amp(i,j)= ( 0.5*( 6.*(M5bar_amp*r_out/J_hinge).^2 + ((F1bar_amp*array_radius)*r_out/I_hinge).^2 ) )^.5;
                % yield stress in steel: 370 MPa-> allow up to 2e8?
                
                rho_st= 8500; %[Kg/m^3]
                area= pi*(r_out^2 -r_in^2);
                mass_hinge(i,j)= rho_st*area*array_radius;
                r_out_mat(i,j)= r_out;
                r_in_mat(i,j)= r_in;
            end
        end
        
        stress_acceptable= hinge_stress_amp( hinge_stress_amp<2e8);
        mass_acceptable= mass_hinge(hinge_stress_amp<2e8);
        r_out_acceptable= r_out_mat(hinge_stress_amp<2e8);
        r_in_acceptable= r_in_mat(hinge_stress_amp<2e8);
        
        [ min_mass, i_min ]= min(mass_acceptable);
        r_out_min= r_out_acceptable(i_min);
        r_in_min= r_in_acceptable(i_min);
        hinge_to_WEC_mass= min_mass/mTMDe;
        
        % check:
        J_hinge= pi*(r_out_min^4-r_in_min^4)/2; % torsional rigidity
        I_hinge= J_hinge/2;
        hinge_stress_check= ( 0.5*( 6.*(M5bar_amp*r_out_min/J_hinge).^2 + ((F1bar_amp*array_radius)*r_out_min/I_hinge).^2 ) )^.5;
        
        
        %figure; surf(mass_hinge)
        %figure;surf(hinge_stress_amp)
        
        % Force steel to be at least 2 cm thick.
        if (r_out_min - r_in_min)< 0.02
            r_in_min=r_out_min - 0.02;
        end
        
        % Outputs
        Hinge.mass= min_mass;
        Hinge.r_out= r_out_min;
        Hinge.r_in= r_in_min;
        Hinge.L= array_radius;
    end
    
    %% Output RAOs
    if RAOrun==1
        %numDOF= 10+3+5*numTMDe; %6 platform + 4 tower + Omega(11) + qDotw(12) + deltaBeta(13) + xTMD(14) + yTMDe(15) + zTMDe(16) + zCol(17) + Pc(18)
        fOut= wVector'./(2*pi); %#ok<*NASGU>
        surgeOut= abs(TF_terms(1,:))';
        heaveOut= abs(TF_terms(3,:))';
        pitchOut= abs(TF_terms(5,:))'.*180./pi;
        stressOut= abs(TF_Stress_surge)'; %[MPa]
        deltaNOut= abs(TF_Nacelle)';
        OmegaOut= abs(TF_terms(11,:))';
        qDotOut= abs(TF_terms(12,:))';
        betaOut= abs(TF_terms(13,:))';
        
        surgeDC= xDC(1);
        heaveDC= 0;
        pitchDC= x5mean*180/pi;
        stressDC= sigmaMean; %[Pa]
        %deltaNDC
        OmegaDC= OmegaMean; %rad/s
        qDotDC= U;
        betaDC= BetaMean; %[DEG]
        
        %Also obtain forcing RAO
        forceOut= X_terms;
                
    if sum( strcmp(geoMod, 'TMDe') )
        WEC_FWT_relativeHeave= S_V33TMDe; 
        WEC_FWT_relativeSurge= S_V11TMDe;
    else
         WEC_FWT_relativeHeave= 0.* TF_terms(1,:); 
         TF_V33TMDe= 0.* TF_terms(1,:); 
    end

        discrete_f= wVector./2/pi; %[Hz]
        discrete_WEC_FWTrelative= (2.*WEC_FWT_relativeHeave.*dw).^.5;
        discrete_heave= (2.* abs(TF_terms(3,:)).^2.*Su .* dw ).^.5;
        discrete_surge= (2.* abs(TF_terms(1,:)).^2.*Su .* dw ).^.5;
        discrete_pitch= (2.* abs(TF_terms(5,:)).^2.*Su .* dw ).^.5 .*180/pi;
        
        % JK 06/16/2019
        % DISCRETE
        figure(2);
        subplot(4,1,1);
        hold on
        plot(discrete_f, discrete_WEC_FWTrelative,'-o');
        box on
        ylabel('Relative FWT-WEC velocity (m/s)')
        subplot(4,1,2);
        hold on
        plot(discrete_f, discrete_heave, '-o');
        box on
        ylabel('FWT Heave (m)')
        subplot(4,1,3);
        hold on
        plot(discrete_f, discrete_surge, '-o');
        box on
        ylabel('FWT Surge (m)')
        subplot(4,1,4);
        hold on
        plot(discrete_f, discrete_pitch, '-o');
        box on
        ylabel('FWT Pitch (deg)')
        xlabel('Frequency (Hz)');
        
%         % CONTINUOUS tf
%         figure(2);
%         subplot(4,1,1);
%         hold on
%         plot(discrete_f, abs(TF_V33TMDe).^2.*Su,'-o');
%         ylabel('Relative FWT-WEC velocity (m/s/m)')
%         subplot(4,1,2);
%         hold on
%         plot(discrete_f, abs(TF_terms(3,:)).^2.*Su, '-o');
%         ylabel('FWT Heave (m/m)')
%         subplot(4,1,3);
%         hold on
%         plot(discrete_f, abs(TF_terms(1,:)).^2.*Su, '-o');
%         ylabel('FWT Surge (m/m)')
%         subplot(4,1,4);
%         hold on
%         plot(discrete_f, abs(TF_terms(5,:).*180./pi).^2.*Su, '-o');
%         ylabel('FWT Pitch (degrees/s)')
        
        %save('FWTAloneTheory_April29');
        
     blah= 1;
        
%         %WEC variables- may not exist
%         if (length(TF_terms(:,1))>13)
%             surgeWECOut= abs(TF_terms(14,:))';
%             heaveWECOut= abs(TF_terms(16,:))';
%             
%             pressureOut= abs(TF_terms(18,:))';
%             
%             powerOut= pressureOut.^2./kWells;
%             
%             relSurgeOut= abs(TF_terms(14,:) - TF_terms(1,:))';
%         end
%         
        
        %UNCOMMENT FOR FLOATING OWC EXPERIMENT
        %     X1Out= abs(forceOut(1,:)' + forceOut(14,:)'+ forceOut(19,:)'); %Rigid in surge
        %     X3Out= abs(forceOut(3,:) - KJy./Lxlinks(1)^2.*TF_terms(3,:) + KJy./Lxlinks(1)^2.*TF_terms(16,:) - KJy./Lxlinks(1).*abs(TF_terms(5,:)) ... %KJy, LxLinks
        %                               -KJy./Lxlinks(2)^2.*TF_terms(3,:) + KJy./Lxlinks(2)^2.*TF_terms(21,:) - KJy./Lxlinks(1).*abs(TF_terms(5,:)) ...
        %                               + K3.*(TF_terms(16,:) + TF_terms(21,:) - 2.*TF_terms(3,:) )    ) ;
        %     X5Out= abs(forceOut(5,:)'); %<-- will add pitch dynamics in simulation
        %     P1Out= abs(TF_terms(18,:))';
        %     pressureOut= P1Out;
        % %     save('linFOWCtheory_May3');    %ACTUAL
        % %     save('NONlinFOWCtheory_May4'); % ACTUAL
        %
        % %      save('theoreticalNONlinFOWC'); %IDEAL
        % %    save('theoreticalLinFOWC'); %IDEAL
        
        
        
        
        % %    UNCOMMENT FOR RIGID1 (3-OWC) EXPERIMENT
        %     X1Out= abs(forceOut(1,:)' + forceOut(14,:)'+ forceOut(19,:)' + forceOut(24,:)'); %Rigid in surge
        %     X3Out= abs(forceOut(3,:)' + forceOut(16,:)'+ forceOut(21,:)' + forceOut(26,:)'...
        %                        + ACol.*(forceOut(18,:)'+ forceOut(23,:)' + forceOut(28,:)')   ); %Rigid in heave. + Water column
        %     X5Out= abs(forceOut(5,:)' + ...
        %                 + forceOut(16,:)'.*-Lxlinks(1)+ forceOut(21,:)'.*-Lxlinks(2) + forceOut(26,:)'.*-Lxlinks(3) ...
        %           + ACol.*(forceOut(18,:)'+ forceOut(23,:)' + forceOut(28,:)') ); %<-- will add pitch dynamics in simulation
        %     P1Out= abs(TF_terms(18,:))';
        %     pressureOut= P1Out;
        %     save('rigid3array_theory_April30');
        
        
        %       UNCOMMENT FOR RIGID2 (BIG CHAMBER) EXPERIMENT
        %     X1Out= abs(forceOut(1,:)' + forceOut(14,:)'); %Rigid in surge. Just 1 chamber
        %     X3Out= abs(forceOut(3,:)' + forceOut(16,:)'+ ACol.*TF_terms(18,:)' ); %Rigid in heave. + Water column
        %     X5Out= abs(forceOut(5,:)'); %<-- will add pitch dynamics in simulation
        %     P1Out= abs(TF_terms(18,:))';
        %     pressureOut= P1Out;
        %     save('rigidBIGarray_theory_May15');
        
        
        %P2Out= abs(TF_terms(23,:))';
        %P3Out= abs(TF_terms(24,:))';
        
        
        
        %   save('nonFOWCtheory_April29');
        
        
    end
    
end




end



