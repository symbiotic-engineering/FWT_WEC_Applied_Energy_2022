%D is the damage inflicted on structure (suggests fraction of life used in lifetime)
function [D_state, n_state]= calculate_Fatigue(sigma, sigmaDot, sigmaMean, p, seaType, materialParam)%, sigmaWind)
%input st dev of stress amplitude and time gradient of stress amplitude [Pa]
%p is frequency of sea/wind state

%n_Life s total cycles in entire device lifetime. Updates during each iteration of calculate_Fatigue
%sigma is in [Pa]

% global sigma_ult m n_Life Tlife. %REMOVE GLOBAL VAR SO PARFOR WORKS
sigma_ult= materialParam(1);
m= materialParam(2);

if strcmp(seaType,'mono')
    a= sigma*2^.5; %a is value of 0-pk stress amplitude [Pa]
    da= 1; %multiple in n_state
    f= 1; %probability distribution of stress amplitudes
else %stress is Gaussian random process --> can use Rayleigh distribution of peaks
    %values of stress magnitude to consider
    a= linspace(0, 8*sigma, 100); %a is value of stress amplitude (0-pk) [Pa]
    da= a(2)-a(1);
    a(1)= []; %NF= inf
    f= a./sigma.^2.*exp(-a.^2./(2.*sigma.^2));%Rayleigh distribution of stress magnitudes
    %sum(f)*da= 1
end

Tz= 2*pi*sigma/sigmaDot; %zero-upcrossing period [s]
Tlife= 20*365.25*24*3600; %Design lifetime in seconds
Tstate= p*Tlife; %time in this state during lifetime [s]

%Tstate/Tz is number of peaks in this state over lifetime
n_state= f.*da.*Tstate./Tz; %vector of # of peaks at each stress amplitude over device lifetime

Sult= sigma_ult; %9.27e8; %8.4e8; %4.15e8; %Steel ultimate stress [Pa]

% Sult
% sigmaMean
% m
NF= .5*((Sult-sigmaMean)./a).^m; %fatigue lifetime for steel at each stress amplitude and mean

D_state= sum(n_state./NF); %Damage due to this state over entire lifetime

% n_Life= n_Life + Tstate/Tz; %n_Life s total cycles in entire device lifetime. Updates during each iteration of calculate_Fatigue
n_state= Tstate/Tz;

%% CHECK FOR MONOCHROMATIC
% DEL_stress= 2*( D_state/(2*n_Life) )^(1/m)*(sigma_ult-abs(sigmaMean)); %this DEL is the pk-pk stess amplitude [Pa]
end



