function P_loss_perc= HVAC_trans_losses(P_Cap_Farm, Power_Wind_Wave_Farm, Wind_hub, dist_to_shore)
%

% FOLLOWING ANALYSIS OF TODOROVIC THESIS

plotOn= 0; %set to 1 to make plots, set to 0 to not plot

% Calculation for farm with one offshore converter station and one onshore converter station

%%%CODE FOR CALCULATING THE LOSSES IN FUNCTION OF WIND%%%%%
%%%%SPEED, I.E. IN FUNCTION OF W.FARM PRODUCTION%%%%%%%%%%%%
f=50;
% R=44.6*1e-6;
% L=0.31*1e-6;
G=0.09*1e-9;
% C=0.26*1e-9; %% all values are per meter!!!%%%%
% B=2*pi*f*C;
% X=2*pi*f*L;

%%%% data from the NEXAN cable brochure 1000mm2, In=l055A,for 132,220KV,
%%%% for 400 KV is Inom=l323 A, 1200 mm2 three voltage levels
In=1055; %%according to the NEXAN data, 1000mm2 cross-sec.%
R_220KV=48*1e-6;
L_220KV=0.378*1e-6;
C_220KV=0.188*1e-9;

% I_nom=1126;

%%%% definition of vector P_loos_wind=f(wind)%%%%%%
% P_loss_wind_220KV=0;

%%%%%% TRANSFORMER DATA %%%%%%%%%%%%%%%%
%%%% in MW%%%%%%%%%%%%%%%%%%%%%%%
P_loss_noload_300MVA_220KV =125*1e-3;
P_loss_noload_400MVA_220KV =80*1e-3;
P_loss_noload_500MVA_220KV =80*1e-3;

%%%% in MW%%%%%%%%%%%%%
P_loss_load_300MVA_220KV=707*1e-3;
P_loss_load_400MVA_220KV =820*1e-3;
P_loss_load_500MVA_220KV=(820-80)*1e-3*((500/400)^2);

Rfe_300MVA_220KV=(220^2)/P_loss_noload_300MVA_220KV;
Rfe_400MVA_220KV=(220^2)/P_loss_noload_400MVA_220KV;
Rfe_500MVA_220KV=(220^2)/P_loss_noload_500MVA_220KV;

%%% in KArnperes %%%%%%%%%%
In_transf_300MVA_220KV=300/(sqrt(3)*220);
In_transf_400MVA_220KV=400/(sqrt(3)*220);
In_transf_500MVA_220KV=500/(sqrt(3)*220);

Rcu_300MVA_220KV=(P_loss_load_300MVA_220KV-P_loss_noload_300MVA_220KV)/(In_transf_300MVA_220KV^2);
Rcu_400MVA_220KV=(P_loss_load_400MVA_220KV-P_loss_noload_400MVA_220KV)/(In_transf_400MVA_220KV^2);
Rcu_500MVA_220KV=(P_loss_load_500MVA_220KV-P_loss_noload_500MVA_220KV)/(In_transf_500MVA_220KV^2);

% initialize for speed
num_hrs= length(Power_Wind_Wave_Farm);
P_loss_wind_220KV=               zeros(dist_to_shore+1,num_hrs); % Cable losses?
P_loss_onshore_transf_220KV=     zeros(dist_to_shore+1,num_hrs); % transformer
P_loss_offshore_transform_220KV= zeros(dist_to_shore+1,num_hrs); % transformer
P_loss_compensation_220KV=       zeros(dist_to_shore+1,num_hrs); % compensation losses
P_losses_220KV=                  zeros(dist_to_shore+1,num_hrs); % distance-saturated P_loss_wind?
P_losses_total_220KV=            zeros(dist_to_shore+1,num_hrs);
P_loss_wind_220KV_perc=          zeros(dist_to_shore+1,num_hrs);
P_losses_total_220KV_perc=       zeros(dist_to_shore+1,num_hrs);

for km=1:dist_to_shore+1 %%% loop for loss calculation of every km-th length
    for hr=1:num_hrs
        L_0=km;
        %%%%%%%% reactive currents for the three voltage levels%%%%%%%%%%%%%
        I_react_220KV=220*(1e+3)*2*pi*f*C_220KV*L_0*1000/(sqrt(3));
        I_react_half_220KV =I_react_220KV /2;
        
        %%%%%% maximal length for the cables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L_max_220KV =In*2*(sqrt(3))/(220*(1e+3)*2 *pi *f*C_220KV*1000);
        
        %%%%%%%%%%%%%%%%%%how many cables are needed for this power????
        number_of_cables_220KV = 3; % Given based on data sources
        
        %%LOSS CALCULATION OF OFFSHORE TRANSFORMER AND COMPENSATION%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if P_Cap_Farm==400
            P_loss_offshore_tr_220KV=((220^2)/Rfe_400MVA_220KV)+Rcu_400MVA_220KV*((Power_Wind_Wave_Farm(hr)/(sqrt(3)*220))^2);
        elseif P_Cap_Farm==500
            P_loss_offshore_tr_220KV=((220^2)/Rfe_500MVA_220KV)+Rcu_500MVA_220KV*((Power_Wind_Wave_Farm(hr)/(sqrt(3)*220))^2);
        elseif P_Cap_Farm==600
            P_loss_offshore_tr_220KV=((220^2)/Rfe_300MVA_220KV)+Rcu_300MVA_220KV*((Power_Wind_Wave_Farm(hr)/2/(sqrt(3)*220))^2);
            P_loss_offshore_tr_220KV=P_loss_offshore_tr_220KV*2;
        elseif P_Cap_Farm== 700
            P_loss_offshore_tr_1_220KV=((220^2)/Rfe_400MVA_220KV)+Rcu_400MVA_220KV*((Power_Wind_Wave_Farm(hr)*4/7/(sqrt(3)*220))^2);
            P_loss_offshore_tr_2_220KV =( (220^2 )/Rfe_300MVA_220KV)+Rcu_300MVA_220KV*( (Power_Wind_Wave_Farm(hr) *3/7 /( sqrt(3 )*220) )^2);
            P_loss_offshore_tr_220KV=P_loss_offshore_tr_1_220KV +P_loss_offshore_tr_2_220KV;
        elseif P_Cap_Farm==800
            P_loss_offshore_tr_220KV=((220^2)/Rfe_400MVA_220KV)+Rcu_400MVA_220KV*((Power_Wind_Wave_Farm(hr)/2/(sqrt(3)*220))^2);
            P_loss_offshore_tr_220KV=P_loss_offshore_tr_220KV*2;
        elseif P_Cap_Farm==900
            P_loss_offshore_tr_1_220KV=((220^2)/Rfe_400MVA_220KV)+Rcu_400MVA_220KV*((Power_Wind_Wave_Farm(hr)*4/9/(sqrt(3)*220))^2);
            P_loss_offshore_tr_2_220KV=( (220^2)/Rfe_500MVA_220KV)+Rcu_500MVA_220KV*( ( Power_Wind_Wave_Farm(hr) *5/9 / ( sqrt(3)*220) )^2);
            P_loss_offshore_tr_220KV=P_loss_offshore_tr_1_220KV +P_loss_offshore_tr_2_220KV;
        else % P_Cap_Farm=> 1000
            P_loss_offshore_tr_1_220KV =((220^2)/Rfe_400MVA_220KV)+Rcu_400MVA_220KV*((Power_Wind_Wave_Farm(hr)*4/10/(sqrt(3)*220))^2 );
            P_loss_offshore_tr_2_220KV=( (220^2 )/Rfe_300MVA_220KV)+Rcu_300MVA_220KV*( ( Power_Wind_Wave_Farm(hr) *3/10/ ( sqrt(3 )*220) )^2 );
            P_loss_offshore_tr_220KV=P_loss_offshore_tr_1_220KV +(P_loss_offshore_tr_2_220KV*2);
        end
        P_loss_offshore_transf_220KV=P_loss_offshore_tr_220KV;
        %%%%%%%%%% COMPENSATION LOSSES CALCULATIONS %%%%%%%%%%%%%%
        Q_tota1_220KV=(2*pi*f)*C_220KV*L_0*1000*(220^2);%%%% in MWs %%%%%%%%%
        P_loss_offshore_comp_220KV=(0.003*2/3)*number_of_cables_220KV*Q_tota1_220KV/2;%% aproximate of
        %%%% the following is just calculation of losses for the one cable and
        %%% at the end it is multiplied with the number of cables to get losses
        %%%% for all cables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%% 220 KV %%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_input_220KV=Power_Wind_Wave_Farm(hr)-P_loss_offshore_transf_220KV-P_loss_offshore_comp_220KV;
        if number_of_cables_220KV> 1
            P_input_220KV=P_input_220KV /number_of_cables_220KV;
        end
        if P_input_220KV<0
            P_input_220KV=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        gama_220KV=sqrt((R_220KV+1j*2*pi*f*L_220KV)*(G+1j*2*pi*f*C_220KV))*1e+3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I1_220KV=P_input_220KV*1e+6/(sqrt(3)*220*1e+3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Zw for three voltage levels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Zw_220KV=sqrt((R_220KV+1j*2*pi*f*L_220KV)/(G+1j*2*pi*f*C_220KV));
        
%         I_nom=1126;
%         Pl_diel=3*0.6;%%%% these values are depending on the cable characteristics, be carefull!!
%         Pl_nom=167.5;
%         Pl_max_l=Pl_nom-Pl_diel;%%% in Watts/meter %%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        tan_vita=0.0003;%% it is const. value for all voltage levels
        Pl_diel_220KV=(2*pi*f)*C_220KV*((220*1e+3)^2)/3*tan_vita; %% in Watts/meter%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pl_max_l_220KV =3*R_220KV*(In^2);
        
%         P_loss4=0;
%         temp=0;
        temp_amb=15; %%%% temperature of ambient 15 Celsius%%%%%%%%%%%%%%%
        alfa_t=0.00393; %%%taken from the paper%%%%%%%%%
        c_alfa=1-alfa_t*(20-temp_amb );
        temp_max=90; %%%%% max temp. of conductor in celsius %%%%%%%%%%%%%
        delta_temp=temp_max-temp_amb;
        
        P_loss4_220KV=0;
%         temp_220KV=0;
        
        %%%%% seting up the line equation for the same distribution of
        %%%%% charging current at both ends%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_at_ends_220KV =I1_220KV-1j*I_react_half_220KV;
        Const_220KV =(I_at_ends_220KV*( cosh(gama_220KV*L_0) )-conj (I_at_ends_220KV) )/sinh(gama_220KV*L_0 );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        I2_220KV= zeros(1, L_0+1);
        I2_abs_220KV= zeros(1, L_0+1);
        temp_220KV= zeros(1, L_0+1);
        for i=1:(L_0+1)
            %%%%% currents calculation for three voltage levels%%%%%%%%%%%%%%%%
            I2_220KV(i)=I_at_ends_220KV*cosh(gama_220KV*(i-1) )-Const_220KV*( sinh(gama_220KV*(i-1)) );
            I2_abs_220KV(i)=abs(I2_220KV(i));
            
            %%%%temp. correction factor temp is for the 4th approximation only %%
            temp_220KV(i)=c_alfa/(c_alfa+alfa_t*delta_temp*(1-((I2_abs_220KV(i)/In)^2)));
            P_loss4_220KV =P_loss4_220KV +(I2_abs_220KV(i)^2)*temp_220KV(i);
        end
        %%%%%%%%%% 4th approximation taking into accout temp. dependance%%%%
        P_loss4_220KV=P_loss4_220KV*Pl_max_l_220KV*((1/In)^2)/L_0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_loss4_220KV =P_loss4_220KV +(Pl_diel_220KV*3);%% "*3" to get Pdiel for all, three cores (phases)
        P_loss4_abs_220KV=abs(P_loss4_220KV);
        %losses for the hr-th wind speed and for all cables, for km-th length in MWs
        P_loss_wind_220KV(km,hr)=P_loss4_abs_220KV*L_0*1e-3*number_of_cables_220KV;
        %%%%%%%%LOSSES OF ONSHORE COMPENSATORS AND TRANSFORMERS%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_input_onshore_trans_220KV =(P_input_220KV*number_of_cables_220KV)-P_loss_wind_220KV(km,hr);
        if P_input_onshore_trans_220KV<0
            P_input_onshore_trans_220KV =0;
        else
        end
        %%AT ONSHORE THERE WILL BE JUST 132/400 AND 220/400 kV ,500 MVA and 300 MVA TRANSFORMERS%%%
%         Rcu_onshore_500MVA_132KV=(808-257)*1e-3/((500/(sqrt(3)*400))^2);
        Rcu_onshore_300MVA_220KV=(820-80)*1e-3/((300/(sqrt(3)*400))^2);
        if P_Cap_Farm<=600
            P_loss_onshore_transf_220KV(km,hr)=(80*1e-3)+Rcu_onshore_300MVA_220KV*( (P_input_onshore_trans_220KV/2/( sqrt(3) *400) )^2 );
            P_loss_onshore_transf_220KV(km,hr)=P_loss_onshore_transf_220KV(km,hr)*2;
        elseif P_Cap_Farm>600 && P_Cap_Farm<=900
            P_loss_onshore_transf_220KV(km,hr)=(80*1e-3)+Rcu_onshore_300MVA_220KV*( (P_input_onshore_trans_220KV /3/( sqrt(3) *400) )^2 );
            P_loss_onshore_transf_220KV(km,hr)=P_loss_onshore_transf_220KV(km,hr)*3;
        elseif P_Cap_Farm>900
            P_loss_onshore_transf_220KV(km,hr)=(80*1e-3)+Rcu_onshore_300MVA_220KV*((P_input_onshore_trans_220KV/4/(sqrt(3)*400))^2);
            P_loss_onshore_transf_220KV(km,hr)=P_loss_onshore_transf_220KV(km,hr)*4;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_loss_offshore_transform_220KV(km,hr)=P_loss_offshore_transf_220KV;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%assumption that onshore and offshore compensations are the same!!!%%%%%%
        P_loss_compensation_220KV (km,hr)=2*P_loss_offshore_comp_220KV;
        %%% END OF THE LOSSES CALCULATION FOR ONSHORE TRANSF. AND COMP.%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Power_Wind_Wave_Farm(hr)==0 %% P_loss_ wind_perc(km,hr) is the cables losses in percent of Wind F.
            %% production for hr-th wind speed and km-th length P_loss_wind_132KV _perc(km,hr)=100;
            P_loss_wind_220KV_perc(km,hr)=100;
        else
            P_loss_wind_220KV_perc(km,hr)=P_loss_wind_220KV(km,hr)/Power_Wind_Wave_Farm(hr)*100;
        end
        %% for the average losses depend. of length for every wind speed%%%
        if km>L_max_220KV
            P_losses_220KV(km,hr)= P_loss_wind_220KV(281,hr);
        else 
            P_losses_220KV(km,hr)= P_loss_wind_220KV(km,hr);
        end
        
        %%%%%TOTAL LOSSES ALL COMPONENTS%%%%%%%%%%%%%%%%%%%%%%%
        P_losses_total_220KV (km,hr)=P_losses_220KV(km,hr)+P_loss_offshore_transform_220KV(km,hr)+...
            P_loss_onshore_transf_220KV(km,hr)+P_loss_compensation_220KV(km,hr);
        if Power_Wind_Wave_Farm(hr)==0 %% P_losses_total_perc(km,hr) is total losses in percent of Wind F.
            %% production for hr-th wind speed and km-th length
            P_losses_total_220KV_perc(km,hr)=100;
        else
            P_losses_total_220KV_perc(km,hr)=P_losses_total_220KV(km,hr)/Power_Wind_Wave_Farm(hr)*100;
        end
    end %% end of the loop for hr=1:length(wf_abs)
end %% end of the loop for km=l :dist

if plotOn == 1
    figure
    plot(Power_Wind_Wave_Farm);
    hold on
    plot(P_losses_total_220KV_perc(end,:),'-b')
    legend('Farm power', 'Power losses');
%     plot(Wind_hub, P_losses_total_220KV_perc(end,:),'-b')
%     ylabel('Percentage of total power losses [%]')
%     xlabel('Wind speed [m/s]')
%     grid
%     axis([0 dist_to_shore 0 L_0])
end

% P_loss = P_losses_total_220KV(end,:);
P_loss_perc = P_losses_total_220KV_perc(end,:)/100;

end




