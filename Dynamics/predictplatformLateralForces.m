function [wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable, A11Table, A15Table, A22Table, A24Table, A33Table, A44Table, A55Table, A66Table, B11Table, B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table]=...
          predictplatformLateralForces (H, zCM, zDraft, z, dz, D)
%HERE, PHASES ARE POSITIVE WHAN LAG, KEEP JONKMAN/WAMIT NOTATION

g= 9.81; %acceleration due to gravity [m/s^2]
rho= 1025; %sea water density [Kg/m^3]

% H= 320; %200 %water depth [m]
% zCM= -89.9155;%depth of platform CM in SWL coordinate system [m]
% zDraft= -120;%[m]
% z
% dz
% D
% % z= linspace(zDraft, 0, 500); %submerged coordinates
% % dz= z(2)-z(1);
% % D= 9.4.*ones(1,length(z));
% % D(z>-12)= -2.9/8.*(z(z>-12)+4)+6.5;
% % D(z>-4)= 6.5; %Diameter for viscous damping
Vsub= pi.*sum(D.^2)./4*dz;

[wForceTable, X1magTableJ, X3magTable, X5magTableJ, X1phaseTableJ, X3phaseTable, X5phaseTableJ, wTable, A11TableJ, A15TableJ, A22TableJ, A24TableJ, A33Table, A44TableJ, A55TableJ, A66Table, B11TableJ, B15TableJ, B22TableJ, B24TableJ, B33Table, B44TableJ, B55TableJ, B66Table]= getWAMIToutputs;
%wForceTable = wTable

wVector= wTable; %linspace(0, 2, 50); %[rad/s]
k= fsolve( @(kk) -wVector.^2+g.*kk.*tanh(kk.*H), .1.*ones(1,length(wVector)));
Vg= (g.*tanh(k.*H)+g.*k.*H.*(sech(k.*H)).^2)./(2.*(g.*k.*tanh(k.*H)).^.5 );


[zz, kk]= ndgrid(z, k);
[DD, kk]= ndgrid(D, k);
[DD, VVg]= ndgrid(D, Vg);

V2d= pi/4*DD.^2;
a11= pi/4*DD.^2;


TF_A1= 1i.*g.*kk.*cosh(kk.*(zz+H))./cosh(kk.*H);
fx= -rho.*(V2d+a11).*TF_A1; %incremental force along platform




A11Table= rho.*sum(a11(:,1))*dz.*ones(1,length(wVector)); %[Kg]
A55Table=  rho.*sum(a11(:,1).*zz(:,1).^2)*dz.*ones(1,length(wVector)); %[Kg]
A15Table= rho.*sum(a11(:,1).*zz(:,1))*dz.*ones(1,length(wVector));
A33= .57*rho*4/24*pi*(D(1)^3).*ones(1,length(wVector)); %.63*rho*4/24*pi*(D(1)^3)
%approximte as semi-infinite circular cylinder with diameter D
    % 1/2 the mass of a circular dick with diameter D because only 1 side wet. 2.24 Recitation 4 Notes
    %MH p. 147 --> limit a-->0. 
A22Table= A11Table.*ones(1,length(wVector));
A24Table= -A15Table.*ones(1,length(wVector));
A44Table= A55Table.*ones(1,length(wVector));



X1= sum(fx, 1)*dz; %[N]
X5= sum(fx.*zz, 1)*dz;

TF_A3= -g.*k.*sinh(k.*(z(1)+H))./cosh(k.*H);
X3= -(rho*Vsub+A33).*TF_A3; %Missing sign change at low frequencies
    
    
%TMDe- approximate Sphere hydrodynamic damping
% B11= k./(8*rho*g.*Vg).*abs(X1).^2 DOESN'T MATCH WAMIT AS WELL
B33= k./(4*rho*g.*Vg).*abs(X3).^2;
    b11_2D= abs(fx).^2./(2*rho*g.*VVg);
    b15_2D= abs(fx).^2./(2*rho*g.*VVg).*zz;
    b55_2D= abs(fx).^2./(2*rho*g.*VVg).*zz.^2;
B11Table= sum(b11_2D, 1)*dz;
B55Table= sum(b55_2D, 1)*dz; %mostly agrees with WAMIT
B15Table= sum(b15_2D, 1)*dz; %I'm overestimating B15 compared to WAMIT

B22Table= B11Table;
B24Table= -B15Table;
B44Table= B55Table;

%KEEP WAMIT 3 COEFFICIENTS AS REMOVE DEPTH OF SPAR.


X1magTable= abs(X1);
X1phaseTable= -angle(X1);
X5magTable= abs(X5);
X5phaseTable= -angle(X5);


end











