function [dates_LMP, LMP_out]= get_hourly_Electricity_Prices

%addpath('Electricity_Prices') % make sure Price data is on path

reprocess_data= 0; % set to 1 to reprocess the data
% else, load processed data from .mat file

if reprocess_data
    % 2017 Data (Maha updated on 3/29/2019):
    %     filenames= {'20170101_20170201_PRC_LMP_DAM_20190207_18_18_17_v1.csv';
    %                 '20170201_20170301_PRC_LMP_DAM_20190207_18_30_53_v1.csv'
    %                 '20170301_20170401_PRC_LMP_DAM_20190207_18_36_29_v1.csv'};
    %
    filenames= {'20170101_20170201_PRC_LMP_DAM_20190329_15_24_48_v1.csv';
        '20170201_20170301_PRC_LMP_DAM_20190329_15_26_05_v1.csv'
        '20170301_20170401_PRC_LMP_DAM_20190329_15_26_53_v1.csv'
        '20170401_20170501_PRC_LMP_DAM_20190329_15_27_52_v1.csv'
        '20170501_20170601_PRC_LMP_DAM_20190329_15_29_34_v1.csv'
        '20170601_20170701_PRC_LMP_DAM_20190329_15_31_02_v1.csv'
        '20170701_20170801_PRC_LMP_DAM_20190329_15_31_57_v1.csv'
        '20170801_20170901_PRC_LMP_DAM_20190329_15_33_48_v1.csv'
        '20170901_20171001_PRC_LMP_DAM_20190329_15_34_38_v1.csv'
        '20171001_20171101_PRC_LMP_DAM_20190329_15_35_28_v1.csv'
        '20171101_20171201_PRC_LMP_DAM_20190329_15_36_20_v1.csv'
        '20171201_20180101_PRC_LMP_DAM_20190329_15_36_57_v1.csv'};
    
    T= [];
    % For each file in filenames, append the data to table T
    for i= 1:length(filenames)
        filename= filenames{i};
        opts= detectImportOptions(filename);
        T = [T; readtable(filename, opts)];
    end
    
    
    % Remove rows for which the 'XML_DATA_ITEM' is 'LMP_CONG_PRC' OR 'LMP_ENE_PRC
    % We only care about 'LMP_PRC'
    T_keep= [];
    for j= 1:height(T)
        if strcmp( T{j, 11}, 'LMP_PRC')
            T_keep= [T_keep; T(j,:)];
        end
    end
    
    % Sort LMP_PRC by the date, and within each date by the hour
    T_final= sortrows(T_keep, [3,4]);
    
    LMP= table2array( T_final(:, 'MW') );
    dates = table2array( T_final(:,'INTERVALSTARTTIME_GMT'));
    
    %    % Repeat for now until get all data for node EUREKAA_6_N001
    % UNNECESSARY BECAUSE GOT ALL DATA (Maha 03/29/2019)
    %     LMP= repmat(LMP, 5, 1);
    %     LMP_out= LMP(1:365*24+1);
    % Awesome, Maha, you're the best! (Jocie 05/04/2019)
    
    % Make sure only get data for 2017
    ind_start = find(strcmp(dates, '2017-01-01T08:00:00-00:00'));
    ind_end = find(strcmp(dates, '2018-01-01T00:00:00-00:00'));
    
    dates_2017 = cell2mat(dates(ind_start:ind_end));
    dates_2017 = datetime([dates_2017(:,1:10) dates_2017(:,12:19)],'InputFormat','yyyy-MM-ddHH:mm:ss');
    LMP_2017 = LMP(ind_start:ind_end);
    
    dates_LMP = datenum(dates_2017)';
    LMP_out = LMP_2017';
    
    % save data to .mat
    save('processed_2017_data', 'dates_LMP', 'LMP_out')
else
    
    load('processed_2017_data', 'dates_LMP', 'LMP_out')
end


end