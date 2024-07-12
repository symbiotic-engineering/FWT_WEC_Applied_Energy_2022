function Demand= demand_fun(date_S)

plotOn= 0; %set to 1 to show plots. Set to 0 to not make plots

%Demand data from CAISO 2017 for the PGE region 

[date, demand]= HourlyLoad_CAISO_PGE_2017;

Demand= interp1(datenum(date), demand, date_S, 'linear', 'extrap');

if plotOn==1
    %Power over year
    figure;
    % plot(hour, demand)
    plot(datenum(date_S), Demand);
    datetick('x','yyyy-mm-dd')
    title('Demand')
    xlabel('hour of year');
    ylabel('demand (MW)');

%     %Power over January 1, June 1, Aug 1
%     hour1= 1:24;
%     demand1= demand(1:24);
% 
%     hour6= 1:24;
%     demand6= demand(3625:3648);
% 
%     hour8= 1:24;
%     demand8= demand(5089:5112);
% 
%     figure(2);
%     plot(hour1, demand1, hour6, demand6, hour8, demand8)
%     title('New England power demand at specific days')
%     xlabel('hour of day');
%     ylabel('demand (kW)');
%     legend('January 1', 'June 1', 'August 1');
end


end