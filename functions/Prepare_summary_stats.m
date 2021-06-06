%% PURPOSE: Run some summary stats for various models
clear all;
addpath('functions')
addpath('data')

%% Prelims
VAR_size = 'HUGE';
date_start = '7/1/1987';  % For full sample results, pu here '7/1/1987'
date_end   = '12/1/2014'; % For full sample results, pu here '12/1/2014'
h = 12;
p = 13; 

series_to_eval  = 1:7;

%% Prepare data
[Y,series,dates]=Prepare_data(VAR_size);
[T,M] = size(Y);
T_thres = round(0.5*T);

% Covert datees from strings to Matlab dates
dates = datenum(char(dates),'dd/mm/yyyy');    

T_beg   = find(dates==datenum(date_start));  
T_end   = find(dates==datenum(date_end));

this_start = T_beg - T_thres + 1;
this_end = T_end - h - T_beg + 1;

% Setup output file
f_id = ['Output/Table_stats_',VAR_size,' (',datestr(date_start,'yyyy.mm'),'--',datestr(date_end,'yyyy.mm'),').xlsx'];


%% Load forecasts from BCVAR and competitors
% BENCH
load([pwd,'/Output/BBENCH_',VAR_size,'.mat']);
% DFM and FAVAR models
load([pwd,'/Output/BCOMP_',VAR_size,'.mat']);

%BCTRVAR, no cov compression
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_50_1_3_2_0.mat']);
fore1     = fore2;
msfe1     = msfe2;
msfe1_ALL = msfe2_ALL;
PL1       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2')

%BCTRVAR, cov compression
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_50_1_3_2_1.mat']);
fore2_1     = fore2;
msfe2_1     = msfe2;
msfe2_1_ALL = msfe2_ALL;
PL2_1       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');

%BCTRVAR with TVP and SV, with cov compression
load([pwd,'/Output/BCTRVAR_TVP_',VAR_size,'_1_50_1_3_2_1.mat']);
fore2_3     = fore2;
msfe2_3     = msfe2;
msfe2_3_ALL = msfe2_ALL;
PL2_3       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');
% BVAR with MINNESOTA priors
load([pwd,'/Output/BVARMINN_',VAR_size,'.mat']);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PL0(PL0<1e-2)     = 1e-2;
% PL1(PL1<1e-2)     = 1e-2;
% PL2_1(PL2_1<1e-2) = 1e-2;
% PL2_3(PL2_3<1e-2) = 1e-2;
% PL3(PL3<1e-2)     = 1e-2;
% PL5(PL5<1e-2)     = 1e-2;
% PL6(PL6<1e-2)     = 1e-2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MSFE ratios for selected series
MSFE_out = [{'Series'},{'Model'},{'Alternative model'},cellstr([repmat('h=',h,1),num2str((1:1:h)')])'];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'DFM'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe4(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'FAVAR1'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe5(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'FAVAR2'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe6(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BVAR MINN'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe3(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BCTRVAR (no cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe1(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BCTRVAR (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe2_1(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BCTRVAR TVP-SV (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe2_3(this_start:this_end,:,:)))'./squeeze(mean(msfe0(this_start:this_end,:,:)))')];

%% WMSFE ratios over selected series
for ii = 1:h
    W = diag(1./diag(cov(Y(T_beg+ii:T_end-h+ii,series_to_eval))));
    tmp(1,ii) = trace(squeeze(sqrt(msfe4(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe4(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
    tmp(2,ii) = trace(squeeze(sqrt(msfe5(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe5(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
    tmp(3,ii) = trace(squeeze(sqrt(msfe6(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe6(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
    tmp(4,ii) = trace(squeeze(sqrt(msfe3(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe3(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
    tmp(5,ii) = trace(squeeze(sqrt(msfe1(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe1(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
    tmp(6,ii) = trace(squeeze(sqrt(msfe2_1(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe2_1(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
    tmp(7,ii) = trace(squeeze(sqrt(msfe2_3(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe2_3(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0(this_start:this_end,ii,:)))');
end
MSFE_out = [MSFE_out;cell(2,size(MSFE_out,2))];
MSFE_out = [MSFE_out;repmat({'WTMSFE (7 series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCTRVAR (no cov)'};{'BCTRVAR (cov)'};{'BCTRVAR TVP-SV (cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];

%% WMSFE ratios over all series
tmp = [];
for ii = 1:h
    W = diag(1./diag(cov(Y(T_beg+ii:T_end-h+ii,:))));
    tmp(1,ii) = trace(squeeze(sqrt(msfe4_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe4_ALL(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))');
    tmp(4,ii) = trace(squeeze(sqrt(msfe3_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe3_ALL(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))');
    tmp(5,ii) = trace(squeeze(sqrt(msfe1_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe1_ALL(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))');
    tmp(6,ii) = trace(squeeze(sqrt(msfe2_1_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe2_1_ALL(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))');
    tmp(7,ii) = trace(squeeze(sqrt(msfe2_3_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe2_3_ALL(this_start:this_end,ii,:)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,:)))');
end
MSFE_out = [MSFE_out;cell(2,size(MSFE_out,2))];
MSFE_out = [MSFE_out;repmat({'WTMSFE (all series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCTRVAR (no cov)'};{'BCTRVAR (cov)'};{'BCTRVAR TVP-SV (cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];

% WMSFE ratios over all series (including financial series)
if isequal(VAR_size,'JURADO')
    tmp = [];
    for ii = 1:h
        W = diag(1./diag(cov(Y(T_beg+ii:T_end-h+ii,1:129))));
        tmp(1,ii) = trace(squeeze(sqrt(msfe4_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe4_ALL(this_start:this_end,ii,1:129)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))');
        tmp(4,ii) = trace(squeeze(sqrt(msfe3_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe3_ALL(this_start:this_end,ii,1:129)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))');
        tmp(5,ii) = trace(squeeze(sqrt(msfe1_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe1_ALL(this_start:this_end,ii,1:129)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))');
        tmp(6,ii) = trace(squeeze(sqrt(msfe2_1_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe2_1_ALL(this_start:this_end,ii,1:129)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))');
        tmp(7,ii) = trace(squeeze(sqrt(msfe2_3_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe2_3_ALL(this_start:this_end,ii,1:129)))')/trace(squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))*W*squeeze(sqrt(msfe0_ALL(this_start:this_end,ii,1:129)))');
    end
    MSFE_out = [MSFE_out;cell(2,size(MSFE_out,2))];
    MSFE_out = [MSFE_out;repmat({'WTMSFE (only FREDMD series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCTRVAR (no cov)'};{'BCTRVAR (cov)'};{'BCTRVAR TVP-SV (cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];
end

%% Diebold and Mariano tests of forecast errors for individual selected series
for ii=1:h
    for jj=series_to_eval
%         pval4_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe4(this_start:this_end,ii,jj),ii,0);
%         pval5_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe5(this_start:this_end,ii,jj),ii,0);
%         pval6_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe6(this_start:this_end,ii,jj),ii,0);
%         pval3_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe3(this_start:this_end,ii,jj),ii,0);
%         pval1_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe1(this_start:this_end,ii,jj),ii,0);
%         pval2_1_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe2_1(this_start:this_end,ii,jj),ii,0);
%         pval2_3_tmp(jj,ii) = test_CR(msfe0(this_start:this_end,ii,jj),msfe2_3(this_start:this_end,ii,jj),ii,0);
        
        [~,~,pval4_tmp(jj,ii),~]   = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe4(this_start:this_end,ii,jj),ii);
        [~,~,pval5_tmp(jj,ii),~]   = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe5(this_start:this_end,ii,jj),ii);
        [~,~,pval6_tmp(jj,ii),~]   = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe6(this_start:this_end,ii,jj),ii);
        [~,~,pval3_tmp(jj,ii),~]   = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe3(this_start:this_end,ii,jj),ii);
        [~,~,pval1_tmp(jj,ii),~]   = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe1(this_start:this_end,ii,jj),ii);
        [~,~,pval2_1_tmp(jj,ii),~] = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe2_1(this_start:this_end,ii,jj),ii);
        [~,~,pval2_3_tmp(jj,ii),~] = test_DM_HLN(msfe0(this_start:this_end,ii,jj),msfe2_3(this_start:this_end,ii,jj),ii);
    end
end

DM_out = [{'Series'},{'Model'},{'Alternative model'},cellstr([repmat('h=',h,1),num2str((1:1:h)')])'];
DM_out = [DM_out;series(series_to_eval)',repmat({'DFM'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval4_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'FAVAR1'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval5_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'FAVAR2'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval6_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BVAR MINN'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval3_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (no cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR TVP-SV (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_3_tmp)];

%% Diebold and Mariano tests of forecast errors over joint forecast errors - selected series
for ii=1:h
    W = diag(1./diag(cov(Y(T_beg+ii:T_end-h+ii,series_to_eval))));
    wsfe0(:,ii) = diag((squeeze((sqrt(msfe0(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe0(this_start:this_end,ii,:)))))');
    wsfe4(:,ii) = diag((squeeze((sqrt(msfe4(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe4(this_start:this_end,ii,:)))))');
    wsfe5(:,ii) = diag((squeeze((sqrt(msfe5(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe5(this_start:this_end,ii,:)))))');
    wsfe6(:,ii) = diag((squeeze((sqrt(msfe6(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe6(this_start:this_end,ii,:)))))');
    wsfe3(:,ii) = diag((squeeze((sqrt(msfe3(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe3(this_start:this_end,ii,:)))))');
    wsfe1(:,ii) = diag((squeeze((sqrt(msfe1(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe1(this_start:this_end,ii,:)))))');
    wsfe2_1(:,ii) = diag((squeeze((sqrt(msfe2_1(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe2_1(this_start:this_end,ii,:)))))');
    wsfe2_3(:,ii) = diag((squeeze((sqrt(msfe2_3(this_start:this_end,ii,:)))))*W*(squeeze((sqrt(msfe2_3(this_start:this_end,ii,:)))))');
end

pval4_tmp = [];
pval5_tmp = [];
pval6_tmp = [];
pval3_tmp = [];
pval1_tmp = [];
pval2_1_tmp = [];
pval2_3_tmp = [];

for ii=1:h
        [~,~,pval4_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe4(this_start:this_end,ii)),ii);
        [~,~,pval5_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe5(this_start:this_end,ii)),ii);
        [~,~,pval6_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe6(this_start:this_end,ii)),ii);
        [~,~,pval3_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe3(this_start:this_end,ii)),ii);
        [~,~,pval1_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe1(this_start:this_end,ii)),ii);
        [~,~,pval2_1_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe2_1(this_start:this_end,ii)),ii);
        [~,~,pval2_3_tmp(1,ii),~] = test_DM_HLN((wsfe0(this_start:this_end,ii)),(wsfe2_3(this_start:this_end,ii)),ii);
end

DM_out = [DM_out;cell(2,size(DM_out,2))];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'DFM'},1,1),repmat({'AR(1)'},1,1),num2cell(pval4_tmp)]];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'FAVAR1'},1,1),repmat({'AR(1)'},1,1),num2cell(pval5_tmp)]];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'FAVAR2'},1,1),repmat({'AR(1)'},1,1),num2cell(pval6_tmp)]];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'BVAR MINN'},1,1),repmat({'AR(1)'},1,1),num2cell(pval3_tmp)]];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'BCTRVAR (no cov)'},1,1),repmat({'AR(1)'},1,1),num2cell(pval1_tmp)]];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'BCTRVAR (cov)'},1,1),repmat({'AR(1)'},1,1),num2cell(pval2_1_tmp)]];
DM_out = [DM_out;[{'WTMSFE (7 series)'},repmat({'BCTRVAR TVP-SV (cov)'},1,1),repmat({'AR(1)'},1,1),num2cell(pval2_3_tmp)]];


copyfile([pwd,'/Input/Table_MSFE_template_new.xlsx'],f_id);
xlswrite(f_id,MSFE_out,'MSFE');
xlswrite(f_id,DM_out,'DM tests (MSFE)');

%% LS differentials for selected series
PL_bench = squeeze(mean(log(PL0)))';

PL_out = [{'Series'},{'Model'},{'Alternative model'},cellstr([repmat('h=',h,1),num2str((1:1:h)')])'];
PL_out = [PL_out;series(series_to_eval)',repmat({'DFM'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL4)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'FAVAR1'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL5)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'FAVAR2'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL6)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'BVAR MINN'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL3)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'BCTRVAR (no cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL1)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'BCTRVAR (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL2_1)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'BCTRVAR TVP-SV (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL2_3)))' - PL_bench)];

%% Joint predictive likelihood for the series to eval across time and forecast horizon

warning('off','MATLAB:nearlySingularMatrix');
disp('Now computing MV PLs...');
for t=1:size(fore0,1)
    for ii=1:h
        yi         = Y(T_thres+t+ii-1,series_to_eval);
        
        this_fore0 = squeeze(fore0(t,:,ii,:));
        this_fore1 = squeeze(fore1(t,:,ii,:));
        this_fore3 = squeeze(fore3(t,:,ii,:));
        this_fore4 = squeeze(fore4(t,:,ii,:));
        this_fore5 = squeeze(fore5(t,:,ii,:));
        this_fore6 = squeeze(fore6(t,:,ii,:));
        this_fore2_1 = squeeze(fore2_1(t,:,ii,:));
        this_fore2_3 = squeeze(fore2_3(t,:,ii,:));
        
        PL0_all(t,ii) = mvnpdf_DP(yi,mean(this_fore0),cov(this_fore0));
        PL1_all(t,ii) = mvnpdf_DP(yi,mean(this_fore1),cov(this_fore1));
        PL3_all(t,ii) = mvnpdf_DP(yi,mean(this_fore3),cov(this_fore3));
        PL4_all(t,ii) = mvnpdf_DP(yi,mean(this_fore4),cov(this_fore4));
        PL5_all(t,ii) = mvnpdf_DP(yi,mean(this_fore5),cov(this_fore5));
        PL6_all(t,ii) = mvnpdf_DP(yi,mean(this_fore6),cov(this_fore6));
        PL2_1_all(t,ii) = mvnpdf_DP(yi,mean(this_fore2_1),cov(this_fore2_1));
        PL2_3_all(t,ii) = mvnpdf_DP(yi,mean(this_fore2_3),cov(this_fore2_3));
    end
end
disp('Done');
warning('on','MATLAB:nearlySingularMatrix');
save([pwd,'/Output/',sprintf('%s_%s.mat','MV-PLS',VAR_size)],'PL0_all','PL1_all','PL3_all','PL4_all','PL5_all','PL6_all','PL2_1_all','PL2_3_all','-mat');

tmp = [];
tmp(1,:) = mean(log(PL4_all)) - mean(log(PL0_all));
tmp(2,:) = mean(log(PL5_all)) - mean(log(PL0_all));
tmp(3,:) = mean(log(PL6_all)) - mean(log(PL0_all));
tmp(4,:) = mean(log(PL3_all)) - mean(log(PL0_all));
tmp(5,:) = mean(log(PL1_all)) - mean(log(PL0_all));
tmp(6,:) = mean(log(PL2_1_all)) - mean(log(PL0_all));
tmp(7,:) = mean(log(PL2_3_all)) - mean(log(PL0_all));

PL_out = [PL_out;cell(2,size(PL_out,2))];
PL_out = [PL_out;repmat({'MVLSD (7 series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCTRVAR (no cov)'};{'BCTRVAR (cov)'};{'BCTRVAR TVP-SV (cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];

%% Diebold and Mariano tests of PL for selected series
pval4_tmp = [];
pval5_tmp = [];
pval6_tmp = [];
pval3_tmp = [];
pval1_tmp = [];
pval2_1_tmp = [];
pval2_3_tmp = [];

for ii=1:h
    for jj=series_to_eval
        [~,~,pval4_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL4(this_start:this_end,ii,jj)),ii);
        [~,~,pval5_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL5(this_start:this_end,ii,jj)),ii);
        [~,~,pval6_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL6(this_start:this_end,ii,jj)),ii);
        [~,~,pval3_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL3(this_start:this_end,ii,jj)),ii);
        [~,~,pval1_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL1(this_start:this_end,ii,jj)),ii);
        [~,~,pval2_1_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL2_1(this_start:this_end,ii,jj)),ii);
        [~,~,pval2_3_tmp(jj,ii),~] = test_DM_HLN(-log(PL0(this_start:this_end,ii,jj)),-log(PL2_3(this_start:this_end,ii,jj)),ii);
    end
end

DM_out = [{'Series'},{'Model'},{'Alternative model'},cellstr([repmat('h=',h,1),num2str((1:1:h)')])'];
DM_out = [DM_out;series(series_to_eval)',repmat({'DFM'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval4_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'FAVAR1'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval5_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'FAVAR2'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval6_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BVAR MINN'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval3_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (no cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR TVP-SV (cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_3_tmp)];

%% Diebold and Mariano tests of PL jointly over selected series
pval4_tmp = [];
pval5_tmp = [];
pval6_tmp = [];
pval3_tmp = [];
pval1_tmp = [];
pval2_1_tmp = [];
pval2_3_tmp = [];

for ii=1:h
        [~,~,pval4_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL4_all(this_start:this_end,ii)),ii);
        %pval4_tmp(1,ii) = NaN;
        [~,~,pval5_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL5_all(this_start:this_end,ii)),ii);
        [~,~,pval6_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL6_all(this_start:this_end,ii)),ii);
        [~,~,pval3_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL3_all(this_start:this_end,ii)),ii);
        [~,~,pval1_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL1_all(this_start:this_end,ii)),ii);
        [~,~,pval2_1_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL2_1_all(this_start:this_end,ii)),ii);
        [~,~,pval2_3_tmp(1,ii),~] = test_DM_HLN(-log(PL0_all(this_start:this_end,ii)),-log(PL2_3_all(this_start:this_end,ii)),ii);
end

DM_out = [DM_out;cell(2,size(DM_out,2))];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'DFM'},1,1),repmat({'AR(1)'},1,1),num2cell(pval4_tmp)]];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'FAVAR1'},1,1),repmat({'AR(1)'},1,1),num2cell(pval5_tmp)]];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'FAVAR2'},1,1),repmat({'AR(1)'},1,1),num2cell(pval6_tmp)]];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'BVAR MINN'},1,1),repmat({'AR(1)'},1,1),num2cell(pval3_tmp)]];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'BCTRVAR (no cov)'},1,1),repmat({'AR(1)'},1,1),num2cell(pval1_tmp)]];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'BCTRVAR (cov)'},1,1),repmat({'AR(1)'},1,1),num2cell(pval2_1_tmp)]];
DM_out = [DM_out;[{'MVLSD (7 series)'},repmat({'BCTRVAR TVP-SV (cov)'},1,1),repmat({'AR(1)'},1,1),num2cell(pval2_3_tmp)]];

xlswrite(f_id,PL_out,'PL');
xlswrite(f_id,DM_out,'DM tests (PL)');
