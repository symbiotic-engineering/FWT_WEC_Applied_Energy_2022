function [Su, wVector, kVector, VgVector, dw]= declareSeaConditions(name, w_p, Hs, H)
%allowable sea inputs: 'Bret', 'Jon', 'white', 'mono'

% H= 320;

wVector= linspace(.01, 1.5, 50);%400); %linspace(.001, 3, 180); %[rad/s]
%wVector= linspace(.01, .15*2*pi, 400); %linspace(.01, .45*2*pi, 200);%400);
dw= wVector(2)-wVector(1);


% BRETSCHNEIDER spectrum
if strcmp(name,'Bret')
    Su= 1.25./4.*(w_p.^4./wVector.^5).*Hs.^2.*exp(-1.25.*(w_p./wVector).^4);  
% JONSWAP spectrum
elseif strcmp(name,'Jon')  
    % From Jonkman thesis PDF p. 47
    Tp= 2*pi/w_p;
    ss= .07.*ones(1,length(wVector));
    ss(wVector>w_p)= .09; 
    func= Tp/Hs^.5;
    if func < 3.6
        gamma= 5;
    elseif func>5
        gamma= 1;
    else
        gamma= exp(5.75-1.15*func);
    end
    Su= 1/(2*pi)*5/16*Hs^2*Tp.*(wVector.*Tp./(2*pi)).^(-5).*exp(-5/4.*(wVector.*Tp./(2*pi)).^(-4)).*(1-.287.*log(gamma)).*gamma.^exp(-.5.*((wVector.*Tp./(2*pi)-1)./ss).^2);

% WHITE NOISE
elseif strcmp(name,'white')  
    %Su= .5*Hs.*ones(1,length(wVector)); %a.^2./(2.*dw).*ones(1,length(wVector));
    Su= (.5*Hs)^2./2./dw.*ones(1,length(wVector));
    
else %strcmp('name','mono')  
    % %Regular wave
    dw= 0.1;
    Su= (Hs/2)^2/2/dw; %Su= a^2/(2*dw)
    wVector= w_p;
    a= (2.*Su.*dw).^.5;

end

g= 9.81; %[m/s^2]

kVector= fsolve( @(kk) -wVector.^2+g.*kk.*tanh(kk.*H), .1.*ones(1,length(wVector)));
VgVector= (g.*tanh(kVector.*H)+g.*kVector.*H.*(sech(kVector.*H)).^2)./(2.*(g.*kVector.*tanh(kVector.*H)).^.5 );

%a= (2.*Su.*dw).^.5; %wave amplitude at each frequency

% % %Remove elements where wave amplitude is negligible
% wVector(a==0)= [];
% kVector(a==0)= [];
% Su(a==0)= [];
% a(a==0)= [];




end







