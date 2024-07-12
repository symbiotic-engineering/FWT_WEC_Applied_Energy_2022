function [wForceTable, X1magTable, X3magTable, X5magTable, X1phaseTable, X3phaseTable, X5phaseTable, wTable A11Table, A15Table, A22Table, A24Table, A33Table, A44Table, A55Table, A66Table, B11Table B15Table, B22Table, B24Table, B33Table, B44Table, B55Table, B66Table]...
    = getWAMIToutputs
%Get hydrodynamic coefficients and forces from WAMIT output files


% coeffs= xlsread('hydroCoeffs.xlsx',1, 'A2:E1021'); 
% 
% coeffs(11:20,:)= []; %remove infinite frequency components
% rho= 1025;
% g= 9.8066;
% 
% T= coeffs(:,1);
% w= 2*pi./T;
% w(1:10)= 0;
% i= coeffs(:,2);
% j= coeffs(:,3);
% A= coeffs(:,4).*rho; %Abar= A/rho
% B= coeffs(:,5).*rho.*w; %Bbar= B/rho/omega
% 
% ii= 1;
% jj= 1;
% while ii<length(w);
%     wTable(jj)= w(ii);
%     A11Table(jj)= A(ii);
%     A15Table(jj)= A(ii+1);
%     A22Table(jj)= A(ii+2);
%     A24Table(jj)= A(ii+3);
%     A33Table(jj)= A(ii+4);
%     A42Table(jj)= A(ii+5);
%     A44Table(jj)= A(ii+6);
%     A51Table(jj)= A(ii+7);
%     A55Table(jj)= A(ii+8);
%     A66Table(jj)= A(ii+9);
%     
%     B11Table(jj)= B(ii);
%     B15Table(jj)= B(ii+1);
%     B22Table(jj)= B(ii+2);
%     B24Table(jj)= B(ii+3);
%     B33Table(jj)= B(ii+4);
%     B42Table(jj)= B(ii+5);
%     B44Table(jj)= B(ii+6);
%     B51Table(jj)= B(ii+7);
%     B55Table(jj)= B(ii+8);
%     B66Table(jj)= B(ii+9);
%     
%     ii= ii+10;
%     jj= jj+1;
% end
% 
% %For now, consider only head-on waves ()
% force= xlsread('waveForce.xlsx',1, 'A2:E22201');
% T= force(:,1);
% w= 2*pi./T;
% angle= force(:,2);
% % i= force(:,3);
% X= force(:,4).*rho.*g;
% phase= force(:,5).*pi./180;
% indShift= 0;%217; %Excel row 218
% 
% ii= 1;
% jj= 1;
% for ii=1:6:length(force)-indShift-5
%     if angle(ii)== -180
%         %ii
%         wForceTable(jj)= w(ii+indShift);
%         X1magTable(jj)= X(ii+indShift);
%         X2magTable(jj)= X(ii+indShift+1);
%         X3magTable(jj)= X(ii+indShift+2);
%         X4magTable(jj)= X(ii+indShift+3);
%         X5magTable(jj)= X(ii+indShift+4);
%         X6magTable(jj)= X(ii+indShift+5);
% 
%         X1phaseTable(jj)= phase(ii+indShift);%+pi; %phase in radians
%         X2phaseTable(jj)= phase(ii+indShift+1);
%         X3phaseTable(jj)= phase(ii+indShift+2);
%         X4phaseTable(jj)= phase(ii+indShift+3);
%         X5phaseTable(jj)= phase(ii+indShift+4);%-pi;
%         X6phaseTable(jj)= phase(ii+indShift+5);
% 
%         jj= jj+1;
%     end
% end
% 
% %add w= 0
%     wForceTable= [0 wForceTable];
%     X1magTable= [X1magTable(1) X1magTable];
%     X3magTable= [X3magTable(1) X3magTable];
%     X5magTable= [X5magTable(1) X5magTable];
%     X1phaseTable= [X1phaseTable(1) X1phaseTable];
%     X3phaseTable= [X3phaseTable(1) X3phaseTable];
%     X5phaseTable= [X5phaseTable(1) X5phaseTable];

load formattedWAMIT_OC3.mat

%Statoil's additional damping for OC3 (Jonkman)
B11Table= B11Table+ 1e5; %[Ns/m]
B22Table= B22Table+ 1e5; 
B33Table= B33Table+ 1.3e5;
B66Table= B66Table+ 13e6;
end
































