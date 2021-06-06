%% Print all tables and figures for February 2016 updated draft

% This now includes results for the BCTRVAR with TVP-SV
% All commented files have not yet been updated to include the new models

%% Appendix tables on transformations
Print_table_transformations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Print_table_MSFE_ratios('MEDIUM',[1 2 3 6 9 12]);
Print_table_MSFE_ratios('LARGE',[1 2 3 6 9 12]);
Print_table_MSFE_ratios('HUGE',[1 2 3 6 9 12]);

% Print_table_MSFE_ratios('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Print_table_MSFE_ratios('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Print_table_MSFE_ratios('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');

%% PL tables
Print_table_PL_diffs('MEDIUM',[1 2 3 6 9 12]);
Print_table_PL_diffs('LARGE',[1 2 3 6 9 12]);
Print_table_PL_diffs('HUGE',[1 2 3 6 9 12]);

% Print_table_PL_diffs('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Print_table_PL_diffs('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Print_table_PL_diffs('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');

%% Relative WTMSFE and MV-PL differential table
Print_table_WTMSFE_ratios([{'MEDIUM'},{'LARGE'},{'HUGE'}],[1 2 3 6 9 12]);
% Print_table_WTMSFE_ratios([{'MEDIUM'},{'LARGE'},{'HUGE'}],[1 2 3 6 9 12],'7/1/1987','12/1/2007');

% %% WCUM SSE Plots
Plot_WCUM_SSE('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_SSE('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_SSE('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');

Plot_WCUM_SSE_zoomed('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_SSE_zoomed('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_SSE_zoomed('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
% 
% Plot_WCUM_SSE('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Plot_WCUM_SSE('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Plot_WCUM_SSE('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');

%% BCVAR TVP-SV tables
Print_table_TVPSV('MEDIUM',[1 2 3 6 9 12]);
Print_table_TVPSV('LARGE',[1 2 3 6 9 12]);
Print_table_TVPSV('HUGE',[1 2 3 6 9 12]);

%% WCUM PL plots
Plot_WCUM_PL('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_PL('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_PL('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');

Plot_WCUM_PL_zoomed('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_PL_zoomed('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');
Plot_WCUM_PL_zoomed('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2014');

% Plot_WCUM_PL('MEDIUM',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Plot_WCUM_PL('LARGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');
% Plot_WCUM_PL('HUGE',[1 2 3 6 9 12],'7/1/1987','12/1/2007');

%% BCVAR TVP-SV plots
%Plot_CUM_SSE_PL_TVPSV([{'MEDIUM'};{'LARGE'};{'HUGE'}],'7/1/1987','12/1/2014');
% Plot_CUM_SSE_TVPSV([{'MEDIUM'};{'LARGE'};{'HUGE'}],'7/1/1987','12/1/2014');
% Plot_CUM_PL_TVPSV([{'MEDIUM'};{'LARGE'};{'HUGE'}],'7/1/1987','12/1/2014');
Plot_CUM_SSE_TVPSV([{'MEDIUM'};{'HUGE'}],'7/1/1987','12/1/2014');
Plot_CUM_PL_TVPSV([{'MEDIUM'};{'HUGE'}],'7/1/1987','12/1/2014');

Plot_CUM_SSE_TVPSV_zoomed([{'MEDIUM'};{'LARGE'};{'HUGE'}],'7/1/1987','12/1/2014');
Plot_CUM_PL_TVPSV_zoomed([{'MEDIUM'};{'LARGE'};{'HUGE'}],'7/1/1987','12/1/2014');

%% Plot on predictive density volatilities
Plot_pred_dens_volatility('MEDIUM',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional plots (for appendices)
Plot_CUM_SSE('MEDIUM',[1 12],'7/1/1987','12/1/2014');
Plot_CUM_SSE('LARGE',[1 12],'7/1/1987','12/1/2014');
Plot_CUM_SSE('HUGE',[1 12],'7/1/1987','12/1/2014');

Plot_CUM_PL('MEDIUM',[1 12],'7/1/1987','12/1/2014');
Plot_CUM_PL('LARGE',[1 12],'7/1/1987','12/1/2014');
Plot_CUM_PL('HUGE',[1 12],'7/1/1987','12/1/2014');

