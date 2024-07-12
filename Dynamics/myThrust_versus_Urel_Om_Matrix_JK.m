function [qW_Matrix, Omega_Matrix, Ct_Matrix, Cq_Matrix, dCt_dq_Matrix, dCq_dq_Matrix, dCt_dOm_Matrix, dCq_dOm_Matrix, Cp_Matrix]= myThrust_versus_Urel_Om_Matrix_JK(bInst)

% bInst= 0; %pitch angle input. Must be integer between 0 and 24 [deg]

rho_air= 1.225; %[Kg/m^3]
R= 63;
A= pi*R^2; %swept area [m^2]

load AeroDyn_CtCqCp2.mat Beta beta TSR tsr Cp Cq Ct 

% figure;
% surf(beta, tsr, Cp); %Each column of Cp is different beta. Each row is different tsr
% xlabel('beta');
% ylabel('tsr');

qW_Vector= linspace(0, 25, 100); %relative wind speed on nacelle [m/s]
dQ= qW_Vector(2)- qW_Vector(1);
Omega_Vector= linspace(0, 30, 200).*2*pi./60; %[rad/s]
dO= Omega_Vector(2) - Omega_Vector(1);

[Omega_Matrix, qW_Matrix]= ndgrid(Omega_Vector, qW_Vector); %each column is different Omega, each row is different qW

TSR_Matrix= Omega_Matrix.*R./qW_Matrix;

%vectors for beta= inputted beta
[~, betaInd]= min(abs((beta-bInst)));

Cq_ref= Cq(:,betaInd);
Ct_ref= Ct(:,betaInd);
Cp_ref= Cp(:,betaInd);

Ct_Matrix= interp1(tsr, Ct_ref, TSR_Matrix); %2D matrix for varied TSR corresponding to Omega, qW, and const beta
    Ct_Matrix(isnan(Ct_Matrix))= 0;
Cq_Matrix= interp1(tsr, Cq_ref, TSR_Matrix); %2D matrix for varied TSR corresponding to Omega, qW, and const beta
    Cq_Matrix(isnan(Cq_Matrix))= 0;
Cp_Matrix= TSR_Matrix.*Cq_Matrix;
    Cp_Matrix(isnan(Cp_Matrix))= 0;
Cp_MatrixCheck= interp1(tsr, Cp_ref, TSR_Matrix); %Agree!
    
dCt_dOm_Matrix= diff(Ct_Matrix, 1,1)./dO;
dCt_dq_Matrix= diff(Ct_Matrix, 1, 2)./dQ; %diff(X,N,DIM)
% dCt_dB_Matrix= diff(Ct_Matrix, 1,3)./dB;
dCq_dOm_Matrix= diff(Cq_Matrix, 1,1)./dO;
dCq_dq_Matrix= diff(Cq_Matrix, 1, 2)./dQ; %diff(X,N,DIM)
% dCq_dB_Matrix= diff(Cq_Matrix, 1,3)./dB;

%Add missing elements to end
dCt_dOm_Matrix(200,:)= 0;
dCt_dq_Matrix(:,100)= 0;
% dCt_dB_Matrix(:,:,25)= 0;
dCq_dOm_Matrix(200,:)= 0;
dCq_dq_Matrix(:,100)= 0;
% dCq_dB_Matrix(:,:,25)= 0;    
    

% %Check is Cp gives power expect
% Drotor= 125.88; %rotor diameter from Jonkman [m]
% % Rrotor= Drotor/2;
% Rrotor= 62.94*cosd(-2.5);
% rho_a= 1.225; %[Kg/m^3]
% S= pi*Drotor^2/4;
% Power_Matrix= .5*rho_a*S.*Cp_Matrix.*qW_Matrix.^3./10^6; 

end












