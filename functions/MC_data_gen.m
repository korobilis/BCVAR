%% Prepare DGPs for Monte Carlo analysis
clear all; close all;
addpath(genpath(pwd));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define some of the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
date_start = '1/1/1960';
date_end   = '12/1/2014';

% Specs for VAR models to simulate
T_sim = 500; %number of obs
p = 2; % number of lags

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the FRED-MD data (updated to 2015:Q1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../../Data/freddata.mat');
[T,N] = size(data);

% Recover dates and turn them in Matlab form
year = floor(dates);
year(dates-floor(dates)==0) = year(dates-floor(dates)==0)-1;
month = round(12*(dates-year));

dates = datenum(year,month,ones(length(year),1));

%% DGP that is fully based on data   

% 1. SMALL. This is a small monetary VAR including the number of employees on 
% non-farm payrolls (EMPL), the consumer price index (CPI) and the Federal 
% Funds Rate (FFR).

vnames_small = [{'PAYEMS'};{'CPIAUCSL'};{'FEDFUNDS'}];

% 2. MEDIUM. addition to the key variables in SMALL, this model includes 
% the index of sensitive material prices (COMM PR) and monetary aggregates: 
% non-borrowed reserves (NBORR RES), total reserves (TOT RES) and M2 money
% stock (M2).  Personal Income INCOME), Real Consumption (CONSUM), Industrial 
% Production (IP), Capacity Utilization (CAP UTIL), Unemployment Rate (UNEMPL), 
% Housing Starts (HOUS START), Producer Price Index (PPI), Personal Consumption 
% Expenditures Price Deflator (PCE DEFL), Average Hourly Earnings (HOUR EARN), 
% M1 Monetary Stock (M1), Standard and Poor’s Stock Price Index (S&P); Yields 
% on 10 year US Treasury Bond (TB YIELD) and effective exchange rate (EXR).

% NOTE: FRED-MD does not include the the index of sensitive material prices (COMM PR)
% (which comes from the BEA), also, the trade-weighted exchange rate (which I used 
% in place of the effective exchange rate) is available in FRED only starting in 1973:1
% so I left it out for now

vnames_med   = [vnames_small;{'M1SL'};{'M2SL'};{'TOTRESNS'};{'NONBORRES'};{'RPI'};{'DPCERA3M086SBEA'};...
               {'INDPRO'};{'CAPUTLB00004S'};{'UNRATE'};{'HOUST'};{'PPIFGS'};{'PCEPI'};...
               {'CES0600000008'};{'S&P 500'};{'GS10'}];
           
% 3. LARGE: This specification includes all available series
vnames_large = names';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare t_beg and t-end vars, will be needed later on
t_beg = find(dates==datenum(date_start));
t_end = find(dates==datenum(date_end));

% Prepare X and Y matrices
X_small = [];
for i=1:length(vnames_small)
    indx = strcmp(vnames_small(i),names');
    disp(' ');
    disp(['Including ',char(names(indx)),' in small VAR']);
    X_small  = [X_small,data(:,indx)];
    nan_indx = isnan(data(:,indx));
    if sum(nan_indx)>0
        disp([char(names(indx)),' is missing on following dates:']);
        disp(datestr(dates(nan_indx),'mmm yyyy'));
    else
        disp([char(names(indx)),' has no missing values']);
    end
end

X_med = [];
for i=1:length(vnames_med)
    indx = strcmp(vnames_med(i),names');
    disp(' ');
    disp(['Including ',char(names(indx)),' in medium VAR']);
    X_med = [X_med,data(:,indx)];
    nan_indx = isnan(data(:,indx));
    if sum(nan_indx)>0
        disp([char(names(indx)),' is missing on following dates:']);
        disp(datestr(dates(nan_indx),'mmm yyyy'));
    else
        disp([char(names(indx)),' has no missing values']);
    end

end

indx_drop = [];
for i=1:length(names)
    disp(' ');
    disp(['Including ',char(names(i)),' in large VAR']);
    nan_indx = isnan(data(:,i));
    if sum(nan_indx)>0 && sum(nan_indx)<24
        disp([char(names(i)),' is missing on following dates:']);
        disp(datestr(dates(nan_indx),'mmm yyyy'));
    elseif sum(nan_indx)== 0
        disp([char(names(i)),' has no missing values']);
    else
        warning([char(names(i)),' has more than 2 years of missing data, drop it from the panel']);
        indx_drop = [indx_drop,i];
    end

end
X_large = data; 
X_large(:,indx_drop) = [];
vnames_large(indx_drop) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate VAR(p) for all three cases, and use estimated parameters to
% simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Y_sim_small,A_OLS_small,SIGMA_OLS_small] = simulate_DGP_actual(dates,X_small,T_sim,p,t_beg,t_end);
[Y_sim_med,A_OLS_med,SIGMA_OLS_med]       = simulate_DGP_actual(dates,X_med,T_sim,p,t_beg,t_end);
[Y_sim_large,A_OLS_large,SIGMA_OLS_large] = simulate_DGP_actual(dates,X_large,T_sim,p,t_beg,t_end);


%% DGP that is sparse and has lots of zeros in the VAR coefficients

% Wait on this, need to figure out a few more things first

%% DGP with factor structure

% Use small simulated VAR from data (from above) to first generate a "loadings 
% matrix" (from U(-1,1)) and project the small VAR to obtain medium and large scale VARs
[Y_sim_med_v2,lam_0_med,lam_1_med]       = simulate_DGP_factor(size(SIGMA_OLS_med,1),dates,X_small,SIGMA_OLS_med,T_sim,1,t_beg,t_end);
[Y_sim_large_v2,lam_0_large,lam_1_large] = simulate_DGP_factor(size(SIGMA_OLS_large,1),dates,X_small,SIGMA_OLS_large,T_sim,1,t_beg,t_end);




