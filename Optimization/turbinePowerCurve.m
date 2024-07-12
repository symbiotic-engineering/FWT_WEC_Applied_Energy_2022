function [U, P]= turbinePowerCurve(P_cap, all_losses)
plotOn= 0; %set to 1 to plot figure. Set to 0 to not plot

%For Siemens SWT-3.6 107. Scaled to 6 MW
%Data compiled by iain Staffell in
% http://www.academia.edu/1489838/Wind_Turbine_Power_Curves

data= [ -1  0
        0   0          %{ wind speed [m/s] Power [kW] }
        1   0
        2   0   
        3   0
        4   80
        5   238
        6   474
        7   802
        8   1234
        9   1773
        10  2379
        11  2948
        12  3334
        13  3515
        14  3577
        15  3594
        16  3599
        17  3600
        18  3600
        19  3600
        20  3600
        21  3600
        22  3600
        23  3600
        24  3600
        25  3600
        26  0
        27  0
        40  0];
%Beyond 25 m/s, turbine cuts-out (stops)

U= data(:,1); %[m/s]
P_orig= data(:,2);

% Scale to 6 MW
P_cap_orig= max(P_orig);
scaling_need= P_cap/P_cap_orig;
P_orig= P_orig.*scaling_need;

% account for losses [Fingersh]
% array_losses = 5/100;
% soiling_losses = 3.5/100;
% all_losses = 12.5/100; % either use array losses and soiling losses, or total losses

% JK Aug. 16, 2019: the Fingersh paper just says soiling losses are 3.5% and
% array losses are 5.0%. Don't see were 12.5% came from. I think we should
% just use 3.5% + 5.0% ~ 8.5% for all_losses.
P = P_orig.*(1-all_losses);

if plotOn==1
    figure(2);
    plot(U,P_orig, U, P);
    xlabel('Wind speed at hub (m/s)');
    ylabel('Turbine output power (kW)');
    title('Wind turbine power curve');
end




end