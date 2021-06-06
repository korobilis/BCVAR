function Compute_forc_stats(DGP_type,VAR_size,p,forc_h,date_est_start,date_for_start,date_for_end)
%% PURPOSE: Load results from MC runs and compute various forecast accuracy metrics

addpath('functions')

%% Get the platform (Unix vs Windows)
OS_tmp = system_dependent('getos');
if ~isempty(strfind(OS_tmp,'Windows'))
    OS = 'Windows';
elseif ~isempty(strfind(OS_tmp,'Linux'))
    OS = 'Linux';
else
    error('Problem getting Operating System - check code')
end
if strmatch(OS,'Windows','exact')
    f_delim = '\';
elseif strmatch(OS,'Linux','exact')
    f_delim = '/';
end


%% Output folder
out_dir = [pwd,f_delim,'Output'];

%% Load the data
load([pwd,f_delim,'Input',f_delim,'MC ',VAR_size,' VAR - ',DGP_type,' based.mat']);
eval(['dates = dates_',VAR_size,';']);
eval(['clear(''dates_',VAR_size,''');']);
eval(['Y_sim = Y_sim_',VAR_size,';']);
eval(['clear(''Y_sim_',VAR_size,''');']);

MC_iter = size(Y_sim,3);

%% Setup start and end of estimation and forecasts, plus a few other things
t_e_beg = min(find(dates>=datenum(date_est_start)));
t_f_beg = min(find(dates>=datenum(date_for_start)));
t_f_end = max(find(dates<=datenum(date_for_end)));

models = [{'BC-SUR v.1'};{'BC-SUR v.2'};{'BC-SUR v.3'};{'BC-VAR v.1'};{'BC-VAR v.2'};{'BC-VAR v.3'};{'BVAR-MINN'};{'BFAVAR'}];
vnames = [{'PAYEMS'};{'CPIAUCSL'};{'FEDFUNDS'}];
% vnames   = [{'PAYEMS'};{'CPIAUCSL'};{'FEDFUNDS'};{'M1SL'};{'M2SL'};{'TOTRESNS'};{'NONBORRES'};{'RPI'};{'DPCERA3M086SBEA'};...
%                {'INDPRO'};{'CAPUTLB00004S'};{'UNRATE'};{'HOUST'};{'PPIFGS'};{'PCEPI'};...
%                {'CES0600000008'};{'S&P 500'};{'GS10'}];

%% Loop over MC iterations and compute RMSFE over whole evaluation period
rmsfe_lab = {};
rmsfe_mat = NaN(0,length(models));
for i=1:16
    file_id_out = [out_dir,f_delim,'MC_output (iter=',num2str(i),'), ',VAR_size,' VAR - ',DGP_type,' based.mat'];
    load(file_id_out);
    for h=1:forc_h
        rmsfe_lab = [rmsfe_lab;[repmat(num2cell(i),length(vnames),1),vnames,repmat(num2cell(h),length(vnames),1)]];
        tmp = NaN(length(vnames),8);
        tmp(:,1) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bcvar1(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,2) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bcvar2(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,3) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bcvar3(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,4) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bcvar4(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,5) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bcvar5(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,6) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bcvar6(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,7) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bvarmin(t_f_beg:t_f_end,:,h)).^2));
        tmp(:,8) = sqrt(nanmean((Y(t_f_beg:t_f_end,:) - Y_forc_bvarfac(t_f_beg:t_f_end,:,h)).^2));
        
        rmsfe_mat = [rmsfe_mat;tmp./repmat(tmp(:,7),1,length(models))];
        
    end
end

%% Write output to excel
col_headers = [{'MC_iter'},{'Varname'},{'forc_h'},models'];
out_cell = [col_headers;[rmsfe_lab,num2cell(rmsfe_mat)]];

f_xls_out = [out_dir,f_delim,'MC_output, ',VAR_size,' VAR - ',DGP_type,' based.xlsx'];
xlswrite(f_xls_out,out_cell,'data')