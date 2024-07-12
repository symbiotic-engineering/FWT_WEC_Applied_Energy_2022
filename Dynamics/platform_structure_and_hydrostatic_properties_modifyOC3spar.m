function [zCMPlatform, mPlatform, IplatformCMx, ITot_Swl_x, C55matplat, C33matplat,...
    mSteelPlatform, mConc_Platform, z_Platform, diameter_Platform,...
    towerCM, mTower, z_Tower, diameter_Tower]=...
        platform_structure_and_hydrostatic_properties_modifyOC3spar...
        (Dbot, D_SWL_platform, D_tower_bot, D_tower_top, mEnd,...
        Lbot, Lmid, LtopSub, Ltop, ...
        tower_height, thickness_platform, thickness_tower, mTMD, zTMD )
    % Adjust ballast as necessary.
    % assume that platform top diameter matches towerbottom diameter
    % Diameter, height, steel thickness parameters

%Z's here are wrt keel, 
% "platform" refers to submerged portion of structure
% "tower" refers to structure of portion above water

%% OC3:
% Ltop= 14;
% Lmid= 8;
% Lbot= 108;
% %LTank= 35.75; %Tank bottom is 4 m from bottom
% Dtop= 6.5;
% Dbot= 9.4;
% tst= .0522; %[m]
% rSphere= 0;
% zSphere= 0;
% mSphere= 0;
% mTMD= 0;

%% Fixed constants
rho_st= 8500; %[Kg/m^3]
rho_conc= 2400; %[Kg/m^3]
rho= 1025;
g= 9.81;

%% Tower geometry

z_startTower= Ltop - LtopSub;
z_Tower= z_startTower:.001:tower_height; % C55 has high sensitivity to dz approx.
    dz= z_Tower(2) - z_Tower(1);
    Dslope= (D_tower_top - D_tower_bot) / (tower_height - z_startTower);
diameter_Tower= D_tower_bot + Dslope.*z_Tower;

dm_Tower= pi*diameter_Tower*thickness_tower*dz*rho_st;

mTower= sum(dm_Tower); % OC3: mTower= 2.5e5 Kg, 78 m tall

% mEnd= 5e5; % Guess for 12 MW Haliade-X. %3.5e5; % Nacelle + Rotor
zEnd= tower_height;

towerCM= sum(dm_Tower .* z_Tower) / mTower; %33.4; %CM of tower wrt SWL

%% Platform geometry 

platform_depth= Lbot + Lmid + LtopSub;
z_keel_platform= -platform_depth;
z_start_bottom_diameter= -LtopSub - Lmid;
z_start_swl_diameter= -LtopSub;


z_Platform= (-abs(platform_depth):.01: z_startTower)'; % C55 has high sensitivity to dz approx.
    dz= z_Platform(2) - z_Platform(1);
    Dslope_mid= ( Dbot - D_SWL_platform) / (z_start_bottom_diameter - z_start_swl_diameter );
diameter_Platform= Dbot * ones( length(z_Platform),1 );
  diameter_Platform( z_Platform > z_start_bottom_diameter)=...
 Dbot + Dslope_mid*( z_Platform(z_Platform > z_start_bottom_diameter ) - z_start_bottom_diameter ); 
diameter_Platform( z_Platform > z_start_swl_diameter)= D_SWL_platform;

% Uncomment to check diameter versus depth
% figure;
% plot(diameter_Platform, z_Platform)

area_Platform= pi*diameter_Platform.^2 ./ 4;
dV_Platform= area_Platform .*dz;
dm_Platform= pi*diameter_Platform*thickness_platform*dz*rho_st;


m_bot_plate= pi * (Dbot/2)^2 * thickness_platform * rho_st;

mSteelPlatform= sum(dm_Platform) + m_bot_plate;

zSteelPlatformCM= ( sum(dm_Platform .* z_Platform ) + m_bot_plate * z_Platform(1) ) / mSteelPlatform; %33.4; %CM of tower wrt SWL



% Platform properties
mSteelPlat= mSteelPlatform; %massTop + MassTop_plate + MassMid + MassBot + MassBot_plate; %Expect 1.70e6
K= z_Platform(1); %-LtotPlatSub; %z-coord of keel below SWL [m] -120
Vsub= sum(dV_Platform); %Vbot + Vmid + LtopSub*pi*(Dtop^2)/4; %total submerged platform volume

%% Concrete Ballast tank

dZbotTank_K= 2; % wrt keel. OC3 used 4 m
Z_bot_Ballast= K + dZbotTank_K;

% Satisfy neutral buoyancy. Place as low as possible in platform.

Fbuoy= rho*g*Vsub;
% Force balance: Fbuoy - g*(mBallast + mPlatform + mTower + mEnd):
mBallast= Fbuoy/g - (mSteelPlatform + mTower + mEnd + mTMD); 

% Potential ballast mass based on platform volume:
dmBallast_vs_z= rho_conc.*area_Platform .* dz; % corresponds to z_Platform

% only consider potential ballast mass starting above botTank constraint
dmBallast_vs_z(z_Platform < Z_bot_Ballast )= []; 
    zBallast_temp= z_Platform;
zBallast_temp(z_Platform < Z_bot_Ballast )= []; 

% Cumulative ballast mass from bottom point versus z
mBallast_temp= cumsum(dmBallast_vs_z);

ind_ZtopBallast= find( mBallast_temp >= mBallast, 1 );
 

% Final ballast parameters
Z_top_Ballast= zBallast_temp( ind_ZtopBallast);
zBallast= z_Platform( z_Platform>= Z_bot_Ballast & z_Platform<= Z_top_Ballast );
diameter_Ballast= diameter_Platform( z_Platform>= Z_bot_Ballast & z_Platform<= Z_top_Ballast );

area_Ballast= area_Platform( z_Platform>= Z_bot_Ballast & z_Platform<= Z_top_Ballast );

dmBallast= rho_conc.*area_Ballast .* dz;

    % OC3: zPlatform= zCM= -89.9 m
mBallast= sum( dmBallast ); %OC3: 5.7663e6
ZBallast_CM= sum(dmBallast .* zBallast) / mBallast;

%% Total platform + Ballast + Tower 
mTotPlat= mSteelPlat + mBallast; % OC3: 7.46633e6

%% COB

% Coordinate transform so origin at SWL

%Vsub= VtopSub + Vmid + Vbot + VSphereSub; %Total submerged volume
B= sum( dV_Platform .* z_Platform )/Vsub; % Center of buoyancy. OC3: B = -62. wrt swl
% K= -(Ltop+Lbot+Lmid); %z-coord of keel below SWL [m] -120
mTot= mSteelPlatform + mTower + mEnd + mBallast + mTMD;
G= ( mSteelPlatform*zSteelPlatformCM + mTower*towerCM + mEnd*zEnd + mBallast*ZBallast_CM + mTMD*zTMD)/mTot; % Total system CM wrt SWL. OC3: -77.5.
% G must be below B to be stable.

Itop_area= pi*( diameter_Platform(end) /2)^4/4;
BM= Itop_area/Vsub; %.0109
KB= B-K; %57.9
KG= G-K; %42.5
% GM= BM + KB-KG; %15.4. Incorrect formula?
GM= BM - G + B;
KM= KG+GM;
M= KM+K; % metacenter

%% OUTPUTS
% ORIGINAL zCM= -89.9155;%depth of platform CM in SWL coordinate system [m]
% platform CM wrt SWL. Accounts for bot plate, platform steel, platform concrete
zCMPlatform= ( mSteelPlatform*zSteelPlatformCM + mBallast*ZBallast_CM)/mTot; % Total system CM wrt SWL. OC3: -77.5.
zDraft= K; %[m]

mSteel_Tower_Platform= mSteelPlatform + mTower + mTMD;
mConc_Platform= mBallast;

% thin hollow cylinder + platform bottom + ballast. Don't include tower here
% IplatformCMx= pi*rho_st*h/12*(3*(rout^4-rin^4)+h^2*(rout^2-rin^2)); thin cylinder?

% moment of inertia: sum: r^2 * dm

mPlatform= mSteelPlatform + mConc_Platform;

    radius_squared= z_Platform.^2 + (diameter_Platform./2).^2 ;  % steel circle's from x axis
Ix_Steel_Sides_Platform_swl= sum( dm_Platform .* radius_squared );

% // Axis theorem
Ix_Steel_Platform_Bot_Plate_swl= 0.25 * m_bot_plate * diameter_Platform(1)^2 + m_bot_plate*z_Platform(1)^2;

% Solid cone shape
% Ix for each incremental disk of concrete
dIx_Conc= 0.25 .* dmBallast .* diameter_Ballast.^2 + dmBallast.*zBallast.^2;
Ix_Conc_Platform_swl= sum(dIx_Conc);

	radius_tower_squared= z_Tower.^2 + (diameter_Tower./2).^2 ;  % steel circle's from x axis
Ix_tower_swl= sum( dm_Tower .* radius_tower_squared );


% Account for platform concrete and steel. Parallel axis theorem
IplatformCMx= Ix_Steel_Sides_Platform_swl + Ix_Steel_Platform_Bot_Plate_swl + Ix_Conc_Platform_swl +...
              mPlatform*zCMPlatform^2  ;
ITot_Swl_x=  Ix_Steel_Sides_Platform_swl + Ix_Steel_Platform_Bot_Plate_swl + Ix_Conc_Platform_swl + Ix_tower_swl;
             

% OC3: 4,229,230,000
   
mTot= mSteelPlatform + mConc_Platform + mTower + mEnd + mTMD;

%Before mooring lines
C55matplat=  mTot*g*GM; %EXPECT OC3 C55= 6.2795e+09 -4999180000 = 1.28e9
                           %C55 positive value comes from buoyancy
                           
% Checks- all C55 should converge                       
Fb= rho*g*Vsub;
W= g*mTot;
Kwp= rho*g*Itop_area;
C55matplat_formula2= -W*G + Fb*B + Kwp;
C55matplat_formula3= -W*G + W*B + Kwp;

err= abs ( (C55matplat_formula3 - C55matplat_formula2)/C55matplat_formula3);

if err > 0.01
    error('Hydrostatic stiffness did not converge to within 1% !')
end
if C55matplat < 0
    error('Unstable in pitch!')
end

% C44matplat= C55matplat;
    Aw= pi*diameter_Platform(end)^2/4;
C33matplat= rho*g*Aw;


% rhogVo= 1025*9.81*V;%6.2795e+09 
end



















