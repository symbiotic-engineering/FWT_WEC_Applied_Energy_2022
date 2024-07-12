function [K77, M77, K99, M99, M71, M75, M91, M95, K79, M79,...
          m11_tower, m55_tower, m15_tower, m51_tower,...
          B77T, B17T, B71T, B75T, B11T, B55T, B15T, B99T,...
          u1TMD, u2TMD]=...
 normalizeEigenShape_Jonkman_really_20160803(zTMD, L_platform_raised)
% L_platform_raised= 10;
% zTMD= 0;

%axial coordinate along tower
L= 77.6; %tower length [m]
z= linspace(0, L, 400);
dz= z(2)- z(1);

r_o= 3.25-.0169.*z;
r_i= 3.223-.0168.*z;
A= pi.*(r_o.^2 - r_i.^2);
E= 210e9; %[Pa] steel
rho_st= 8500; %structure density [Kg/m^3]
I= pi./4.*(r_o.^4 - r_i.^4); %[m^4]

% L_SWL_TB= L_platform_raised; %COORD SYSTEM AT SWL  %89.9155+10; %distance from SWL to tower base [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mmTower=249718; %347460; %[Kg]
%     mmPlatform=7466330; %[Kg]
    mmNacelleHub= 56780+240000; %240000; %56780+240000; %hub + nacelle [Kg]
    mmBlades= 3*17740; %rotor mass
    IxxRotor= 2*11776047; %Result of 3 BALANCED blades
    IHub= 115926;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endMass= mmNacelleHub;%+mmRotor;%+mmRotor;%+mmRotor;%COMMENTED JUNE 17, 2016 + mmRotor; % %Nacelle+hub end mass [Kg]
endMassRigid= (mmNacelleHub + mmBlades);
endInertia= IHub+IxxRotor;% + endMass.*(L_SWL_TB+L).^2;% + mmNacelle <-- accounted for in endMass
endInertiaRigid= IxxRotor + IHub;

% %Jonkman Mode 1, Fore-Aft
A1=  0.8689/L^2;%  coefficient of x^2 term
B1=  0.2205/L^3;%  coefficient of x^3 term
C1= -0.0908/L^4;%  coefficient of x^4 term
D1=  0.1167/L^5;%  coefficient of x^5 term
E1=  -0.1154/L^6;% coefficient of x^6 term

% % %Jonkman Mode 2, Fore-Aft
A2=  42.5859/L^2;%  coefficient of x^2 term
B2=  -18.6419/L^3;%  coefficient of x^3 term
C2= -20.3570/L^4;%  coefficient of x^4 term
D2=  -23.2686/L^5;%  coefficient of x^5 term
E2=  20.6816/L^6;% coefficient of x^6 term

%%%%%%%%%%%%%%%%%%%% Integral coefficients %%%%%%%%%%%%%%%%
u1=           A1.*z.^2 + B1.*z.^3    + C1.*z.^4 + D1.*z.^5    +     E1.*z.^6;
u1_grad_1= 2.*A1.*z + 3.*B1.*z.^2 + 4.*C1.*z.^3 + 5.*D1.*z.^4 +  6.*E1.*z.^5;
u1_grad_2= 2.*A1 +    6.*B1.*z +   12.*C1.*z.^2 +20.*D1.*z.^3 + 30.*E1.*z.^4;
% u1_grad_3=            6.*B1 +      24.*C1.*z +   60.*D1.*z.^2 +120.*E1.*z.^3;
% u1_grad_4=                         24.*C1    +  120.*D1.*z    +360.*E1.*z.^2;  

u2=           A2.*z.^2 + B2.*z.^3  +   C2.*z.^4 +    D2.*z.^5 +     E2.*z.^6;
u2_grad_1= 2.*A2.*z + 3.*B2.*z.^2 + 4.*C2.*z.^3 + 5.*D2.*z.^4 +  6.*E2.*z.^5;
u2_grad_2= 2.*A2 +    6.*B2.*z +   12.*C2.*z.^2 +20.*D2.*z.^3 + 30.*E2.*z.^4;
% u2_grad_3=            6.*B2 +      24.*C2.*z +   60.*D2.*z.^2 +120.*E2.*z.^3;
% u2_grad_4=                         24.*C2    +  120.*D2.*z    +360.*E2.*z.^2;  



M71= rho_st*sum(A.*u1.*dz) + endMass*u1(end); %in EOM_1 row. coefficient of xp1''
M75= rho_st.*sum(A.*u1.*(L_platform_raised+z).*dz)     + endMass*u1(end)*(L_platform_raised+z(end)) + endInertia*u1_grad_1(end); %last term NEW
M77= rho_st*sum(A.*u1.^2.*dz) + endMass*u1(end)^2 + endInertia*u1_grad_1(end)^2; %last term NEW
M79= rho_st*sum(A.*u1.*u2.*dz) + endMass*u1(end)*u2(end) + endInertia*u1_grad_1(end)*u2_grad_2(end); %last term NEW

K77= E*sum(I.*u1_grad_2.^2.*dz); %K77
K79= E*sum(I.*u2_grad_2.*u1_grad_2.*dz); %negligible? Best-fit of mode 2 shape off? Large mode 2 magnitude causes large magnitude here?


% M77= .6*M77;
% K77= .6*K77;
M71= M71;%#ok<ASGSL> %1.87*M71;
M75= M75;%#ok<ASGSL> %1.87*M75;
% M71= 1.67*M71;
% M75= 1.67*M75;

% % K77_check= E*sum(I.*u1_grad_4.*u1.*dz); %K77? Error due to not being true eigenshape or numerical error?
% w1= (K77/M77)^.5;
% f1= w1/(2*pi); %w= 2pif. f= w/(2pi). Expect 0.426 Hz
% % f1_check= (K1_check/M77)^.5/(2*pi);


M91= rho_st*sum(A.*u2.*dz) + endMass*u2(end);
M95= rho_st.*sum(A.*u2.*(L_platform_raised+z).*dz) + endMass*u2(end)*(L_platform_raised+z(end)) + endInertia*u2_grad_2(end); %New
M99= rho_st*sum(A.*u2.^2.*dz) + endMass*u2(end)^2 + endInertia*u2_grad_1(end)^2;

K99= E*sum(I.*u2_grad_2.^2.*dz);

% M99= .6*M99;
% K99= .6*K99;
M91= M91; %#ok<ASGSL> %1.87*M91;
M95= M95;%#ok<ASGSL> %1.87*M95;
% M91= 1.67*M91;
% M95= 1.67*M95;

% % K99_check= E*sum(I.*u2_grad_4.*u2.*dz);
% % M99= rho_st*sum(A.*u2.^2.*dz) + endMass*u2(end)^2;
% w2= (K99/M99)^.5;
% f2= w2/(2*pi); %w= 2pif. f= w/(2pi). Expect 3.743 Hz
% % f2_check= (K2_check/M99)^.5/(2*pi);


%Rigid mass properties
m11_tower= rho_st*sum(A)*dz + endMassRigid; %endMass;
m55_tower= rho_st*sum(A.*(L_platform_raised+z).^2)*dz + endMassRigid*(L_platform_raised+L)^2 +  endInertiaRigid;
m15_tower= rho_st*sum(A.*(L_platform_raised+z))*dz + endMassRigid*(L_platform_raised+L);
m51_tower= rho_st*sum(A.*(L_platform_raised+z))*dz + endMassRigid*(L_platform_raised+L);


% % STRESS
% u1factor= 1/u1(end);%multiply stress by this to u1(end)= 1.
% u2factor= 1/u2(end);
% sigma_Root1= E*r_o(1)*( -u1_grad_2(1))/(10^6)*u1factor; %sigma= EI*kappa at root. Coefficient of alpha_1 [MPa]
% sigma_Root2= E*r_o(1)*( -u2_grad_2(1))/(10^6)*u2factor; %sigma= EI*kappa at root. Coefficient of alpha_2 [MPa]

%NONDIAGONAL TERMS (assume numerical errors in gradients? Negligible)
% K79_check= E*sum(I.*u2_grad_4.*u1.*dz); %In EOM_1 row. Coefficient of alpha1
% K97_check= E*sum(I.*u1_grad_4.*u2.*dz);%In EOM_2 row Coefficienct of alpha2



% %FORCING ON TOWER
% M71= rho_st*sum(A.*u1.*dz) + endMass*u1(end); %in EOM_1 row. coefficient of xp1''
% M75= rho_st.*sum(A.*u1.*(L_platform_CM_TB+z).*dz) + endMass*u1(end)*(L_platform_CM_TB+z(end)); %in EOM_1 row. coefficient of xp5''

% M91= rho_st*sum(A.*u2.*dz) + endMass*u2(end); %in EOM_2 row. coefficient of xp1''
% M95= rho_st.*sum(A.*u2.*(L_platform_CM_TB+z).*dz) + endMass*u2(end)*(L_platform_CM_TB+z(end)); %in EOM_1 row. coefficient of xp5''


%Set to 0 for rigid tower
%m17_tower= rho_st*sum(A.*u1)*dz + endMass*u1(end);
% m57_tower= rho_st*sum(A.*u1.*(L_platform_CM_TB+z))*dz + endMass*(L_platform_CM_TB+L)*u1(end);
% 
% m19_tower= rho_st*sum(A.*u2)*dz + endMass*u2(end);
% m59_tower= rho_st*sum(A.*u2.*(L_platform_CM_TB+z))*dz + endMass*(L_platform_CM_TB+L)*u2(end);

 
% dK= E.*I.*dz;
% dm= rho_st.*A.*dz;
% b= .001.*2.*dm.^.5; %(dK.*dm).^.5; %incremental damping

B77T= 0;%.01*(K77*M77)^.5; %sum(b.*u1);
B17T= 0; %B77T;
B71T= 0; %B17T;
% B57T= 0; %sum(b.*u1.*(L_SWL_TB+z));
B75T= 0; %B57T;
B11T= 0; %B77T;
B55T= 0;%L_SWL_TB^2*B11T;
B15T= 0;%L_SWL_TB*B11T;

B99T= 0;%.01*(K99*M99)^.5;



%TMD ASSUME TOWER IS 10 M ABOVE WATERLINE
u1TMD= interp1(z, u1, zTMD-L_platform_raised); %coefficient for alpha1 displacement at TMD height
u2TMD= interp1(z, u2, zTMD-L_platform_raised);



end








