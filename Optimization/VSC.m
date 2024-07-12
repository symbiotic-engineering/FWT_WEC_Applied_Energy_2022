function [P_wf, P_loss, P_in_CS, loss_CS, P_in_cab, P_loss_cab]= VSC(P_Farm)
%Table of [Wind power before offshore substation, power delivered to grid]

plotOn= 0; %set to 1 to make plots, set to 0 to not plot

% Calculation for 500 MW wind farm with one offshore converter station and one onshore converter station

% %Farm C has 500 MW capacity
% [pow_in_C,loss_per_conv_C]=MW500_loss_fun; %Loss curve of a 500 MW converter station
% %[power out, power in, losses as a percent] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% our code for farm generated power versus wind speed
% [wind, P_wf]= turbinePowerCurve; %JK edit May 5, 2015 
% Nturb= 139; %number of turbines for 500 MW farm
% P_wf= P_wf*Nturb./1000; %Farm output power corresponding to wind speed [MW]

P_wf= linspace(0, P_Farm, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for the cable based on Brakelmann
l_c= 80; %Farm distance to shore. [km]
l_cable=3*l_c; %Total cable length [km]
v_rated_cable=3*220; % [kV]
T_amb=15; %ambient [water] temperature [C]
T_op_max=70; %maximum allowable operating temperature for cable [C]

%Cable parameters from Standards- See Barberis PDF p. 85
r0_C= 0.0138; % scalable resistance per length. From Standards with section of 1300 mm2 [Ohm/m]
R_20_C= r0_C*l_cable; % total cable resistance [ohm]
alfa20_C= 3.93e-3; % constant mass temperature coefficient at 20C from Standards [1/K];
c_alfa_C= 1-alfa20_C*(20-T_amb);
Dthetamax= T_op_max-T_amb; % maximum allowable difference between cable temp and ambient temp [C]
c_m_C= 1+alfa20_C*(Dthetamax+T_amb-20);
i_nom_C= 1.677; % nominal current of cable [kA]
Plmax_C= R_20_C*c_m_C*i_nom_C^2; % Maximum power lost in cable when operating at max rating
    

%CS parameters
% JK: I created these variables to use Barberis's estimate of losses in 500 MW CS
% without doing his intermediate computations
%power into CS [MW]

[P_in_CS, per_loss_CS]= estimate_500MW_CS_Losses(P_Farm);

loss_CS= per_loss_CS.*P_in_CS; %power lost in converter station is %power loss * power into station


for ii=1:length(P_wf) %for each wind speed:

    %1) First, power is lost in offshore substation 1
    P_loss_stat1(ii)= interp1(P_in_CS, loss_CS, P_wf(ii));
    %[P_loss_stat1(ii),loss_stat1(ii)]=loss_stat_500(P_wf(ii),pow_in_C,loss_per_conv_C); %losses in station 1

    %2) Losses in cables--> calculate cable reisistance, rr, for this power flow
    P_in_cab(ii)=P_wf(ii)-P_loss_stat1(ii); %Power flow through cable [from CS1 to CS2] [MW]
    i_in_cab(ii)=P_in_cab(ii)/v_rated_cable; %current in cable
        %Cable parameter from Brakelman
    v_theta(ii)=(c_alfa_C)/(c_alfa_C+alfa20_C*Dthetamax*(1-(i_in_cab(ii)/i_nom_C)^2)); %cable parameter
    P_loss_cab(ii)=Plmax_C*(i_in_cab(ii)/i_nom_C)^2*v_theta(ii); %power lost in cable
    r(ii)=P_loss_cab(ii)/i_in_cab(ii)^2; %cable total resistance 
    
    %3) Power losses in onshore susbtation 2
    P_in_stat2(ii)=P_in_cab(ii)-P_loss_cab(ii); %power into station 2
    P_loss_stat2(ii)= interp1(P_in_CS, loss_CS, P_in_stat2(ii));
    

    %Power delivered to grid
    P_out(ii)=P_in_stat2(ii)-P_loss_stat2(ii); %power out of station 2
    
    %total station losses
    P_loss_stat(ii)=P_loss_stat1(ii)+P_loss_stat2(ii); %losses of both stations

    if P_out(ii)<=0, %if losses exceed power supplied by turbines
        loss_per_tot(ii)=100; %power lost is 100%
        P_out(ii)=0; %output power is 0
    else
        loss_per_tot(ii)=(P_wf(ii)-P_out(ii))/P_wf(ii)*100; %percent turbine power lost by transmission system
    end
end

%Output vectors
% loss_per_tot; %percent of wind farm turbine power that is lost
P_loss=loss_per_tot.*P_wf./100; %Magnitude of turbine power lost [MW]

if plotOn== 1
    figure(1)
    subplot(2,1,1)
    plot(P_wf, P_wf-P_out, P_wf, P_loss_stat1, P_wf, P_loss_cab, P_wf, P_loss_stat2);
    title('Power losses of transmission system- NO STORAGE')
    xlabel('Wind turbine power (MW)')
    ylabel('Power output by transmission system (MW)');
    legend('Total losses', 'Offshore susbstation', 'Cable', 'Onshore substation');
    xlim([0 max(P_wf)]);
    subplot(2,1,2);
    plot(P_wf, loss_per_tot, P_wf, [100 P_loss_stat1(2:end)./P_wf(2:end).*100], P_wf, [100  P_loss_cab(2:end)./P_wf(2:end).*100], P_wf, [100 P_loss_stat2(2:end)./P_wf(2:end).*100]);
    title('Percent wind turbine power lost by transmission system- NO STORAGE')
    xlabel('Wind turbine power (MW)')
    ylabel('% Power');
    xlim([0 max(P_wf)]);
end
P_in_cab(1)= [];
P_loss_cab(1)= [];

end




