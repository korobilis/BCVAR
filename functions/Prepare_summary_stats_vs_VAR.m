%% PURPOSE: Run some summary stats for various models
clear all;
addpath('functions')
addpath('data')

%% Prelims
VAR_size = 'LARGE';
h = 12;
p = 13; 

series_to_eval  = 1:7;

%% Prepare data
[Y,series,dates]=Prepare_data(VAR_size);
[T,M] = size(Y);
T_thres = round(0.5*T);

% Covert datees from strings to Matlab dates
dates = datenum(char(dates),'dd/mm/yyyy');    

%% Load forecasts from BCVAR and competitors
% AR(1) benchmark
load([pwd,'/Output/BBENCH_',VAR_size,'.mat']);
fore0_1     = fore0;
msfe0_1     = msfe0;
msfe0_1_ALL = msfe0_ALL;
PL0_1       = PL0;
% Small VAR benchmark
load([pwd,'/Output/BBENCH_VAR_',VAR_size,'.mat']);
% DFM and FAVAR models
load([pwd,'/Output/BCOMP_',VAR_size,'.mat']);
% BCVAR
msfe1     = msfe0;
msfe1_ALL = msfe0_ALL;
PL1       = PL0;
fore1     = fore0;

%load([pwd,'/Output/BCVAR_',VAR_size,'_1_1.mat']);
%BCTRVAR, no cov compression
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_50_1_3_2_0.mat']);
fore2_1     = fore2;
msfe2_1     = msfe2;
msfe2_1_ALL = msfe2_ALL;
PL2_1       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');
%BCTRVAR, with cov compression
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_50_1_3_2_1.mat']);
fore2_3     = fore2;
msfe2_3     = msfe2;
msfe2_3_ALL = msfe2_ALL;
PL2_3       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');

if ~isequal(VAR_size,'JURADO')
    load([pwd,'/Output/BVARMINN_',VAR_size,'.mat']);
else
    msfe3     = msfe0;
    msfe3_ALL = msfe0_ALL;
    PL3       = PL0;
    fore3     = fore0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PL0(PL0<1e-2)     = 1e-2;
PL0_1(PL0_1<1e-2) = 1e-2;
PL1(PL1<1e-2)     = 1e-2;
PL2_1(PL2_1<1e-2) = 1e-2;
PL2_3(PL2_3<1e-2) = 1e-2;
PL3(PL3<1e-2)     = 1e-2;
PL5(PL5<1e-2)     = 1e-2;
PL6(PL6<1e-2)     = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prep LHS series
for irep = T_thres:T-h
    for ii = 1:h
        yi(irep-T_thres+1,:,ii) = Y(irep+ii,series_to_eval);
    end
end

% MSFE ratios for selected series
MSFE_out = [{'Series'},{'Model'},{'Alternative model'},cellstr([repmat('h=',h,1),num2str((1:1:h)')])'];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'DFM'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe4))'./squeeze(mean(msfe0))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'FAVAR1'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe5))'./squeeze(mean(msfe0))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'FAVAR2'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe6))'./squeeze(mean(msfe0))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BVAR MINN'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe3))'./squeeze(mean(msfe0))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BCVAR'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe1))'./squeeze(mean(msfe0))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe2_1))'./squeeze(mean(msfe0))')];
MSFE_out = [MSFE_out;cell(1,size(MSFE_out,2))];
MSFE_out = [MSFE_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3 cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(msfe2_3))'./squeeze(mean(msfe0))')];

% weighted trace mean squared forecast error ratios over the 7 series
for ii = 1:h
    W = diag(1./diag(cov(Y(T_thres+ii:T-h+ii,series_to_eval))));
    tmp(1,ii) = trace(squeeze(sqrt(msfe4(:,ii,:)))*W*squeeze(sqrt(msfe4(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
    tmp(2,ii) = trace(squeeze(sqrt(msfe5(:,ii,:)))*W*squeeze(sqrt(msfe5(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
    tmp(3,ii) = trace(squeeze(sqrt(msfe6(:,ii,:)))*W*squeeze(sqrt(msfe6(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
    tmp(4,ii) = trace(squeeze(sqrt(msfe3(:,ii,:)))*W*squeeze(sqrt(msfe3(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
    tmp(5,ii) = trace(squeeze(sqrt(msfe1(:,ii,:)))*W*squeeze(sqrt(msfe1(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
    tmp(6,ii) = trace(squeeze(sqrt(msfe2_1(:,ii,:)))*W*squeeze(sqrt(msfe2_1(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
    tmp(7,ii) = trace(squeeze(sqrt(msfe2_3(:,ii,:)))*W*squeeze(sqrt(msfe2_3(:,ii,:)))')/trace(squeeze(sqrt(msfe0(:,ii,:)))*W*squeeze(sqrt(msfe0(:,ii,:)))');
end
MSFE_out = [MSFE_out;cell(2,size(MSFE_out,2))];
MSFE_out = [MSFE_out;repmat({'WTMSFE (7 series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCVAR'};{'BCTRVAR (opt 3)'};{'BCTRVAR (opt 3 cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];

% weighted trace mean squared forecast error ratios over all series
tmp = [];
for ii = 1:h
    W = diag(1./diag(cov(Y(T_thres+ii:T-h+ii,:))));
    tmp(1,ii) = trace(squeeze(sqrt(msfe4_ALL(:,ii,:)))*W*squeeze(sqrt(msfe4_ALL(:,ii,:)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,:)))');
    tmp(4,ii) = trace(squeeze(sqrt(msfe3_ALL(:,ii,:)))*W*squeeze(sqrt(msfe3_ALL(:,ii,:)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,:)))');
    tmp(5,ii) = trace(squeeze(sqrt(msfe1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe1_ALL(:,ii,:)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,:)))');
    tmp(6,ii) = trace(squeeze(sqrt(msfe2_1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe2_1_ALL(:,ii,:)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,:)))');
    tmp(7,ii) = trace(squeeze(sqrt(msfe2_3_ALL(:,ii,:)))*W*squeeze(sqrt(msfe2_3_ALL(:,ii,:)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,:)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,:)))');
end
MSFE_out = [MSFE_out;cell(2,size(MSFE_out,2))];
MSFE_out = [MSFE_out;repmat({'WTMSFE (all series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCVAR'};{'BCTRVAR (opt 3)'};{'BCTRVAR (opt 3 cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];

% weighted trace mean squared forecast error ratios over all series
if isequal(VAR_size,'JURADO')
    tmp = [];
    for ii = 1:h
        W = diag(1./diag(cov(Y(T_thres+ii:T-h+ii,1:129))));
        tmp(1,ii) = trace(squeeze(sqrt(msfe4_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe4_ALL(:,ii,1:129)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))');
        tmp(4,ii) = trace(squeeze(sqrt(msfe3_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe3_ALL(:,ii,1:129)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))');
        tmp(5,ii) = trace(squeeze(sqrt(msfe1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe1_ALL(:,ii,1:129)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))');
        tmp(6,ii) = trace(squeeze(sqrt(msfe2_1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe2_1_ALL(:,ii,1:129)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))');
        tmp(7,ii) = trace(squeeze(sqrt(msfe2_3_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe2_3_ALL(:,ii,1:129)))')/trace(squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))*W*squeeze(sqrt(msfe0_1_ALL(:,ii,1:129)))');
    end
    MSFE_out = [MSFE_out;cell(2,size(MSFE_out,2))];
    MSFE_out = [MSFE_out;repmat({'WTMSFE (only FREDMD series)'},7,1),[{'DFM'};{'FAVAR1'};{'FAVAR2'};{'BVAR MINN'};{'BCVAR'};{'BCTRVAR (opt 3)'};{'BCTRVAR (opt 3 cov)'}],repmat({'AR(1)'},7,1),num2cell(tmp)];
end

% Diebold and Mariano tests
for ii=1:h
    for jj=series_to_eval
        pval4_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe4(:,ii,jj),ii,0);
        pval5_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe5(:,ii,jj),ii,0);
        pval6_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe6(:,ii,jj),ii,0);
        pval3_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe3(:,ii,jj),ii,0);
        pval1_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe1(:,ii,jj),ii,0);
        pval2_1_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe2_1(:,ii,jj),ii,0);
        pval2_3_tmp(jj,ii) = test_CR(msfe0(:,ii,jj),msfe2_3(:,ii,jj),ii,0);
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
DM_out = [DM_out;series(series_to_eval)',repmat({'BCVAR'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3 cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_3_tmp)];

f_id = ['Output/Table_MSFE_',VAR_size,'(vs VAR).xlsx'];
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
PL_out = [PL_out;series(series_to_eval)',repmat({'BCVAR'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL1)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL2_1)))' - PL_bench)];
PL_out = [PL_out;cell(1,size(PL_out,2))];
PL_out = [PL_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3 cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(squeeze(mean(log(PL2_3)))' - PL_bench)];

% Diebold and Mariano tests of forecast errors
pval4_tmp = [];
pval5_tmp = [];
pval6_tmp = [];
pval3_tmp = [];
pval1_tmp = [];
pval2_1_tmp = [];
pval2_3_tmp = [];

for ii=1:h
    for jj=series_to_eval
        pval4_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL4(:,ii,jj)),ii,0);
        pval5_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL5(:,ii,jj)),ii,0);
        pval6_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL6(:,ii,jj)),ii,0);
        pval3_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL3(:,ii,jj)),ii,0);
        pval1_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL1(:,ii,jj)),ii,0);
        pval2_1_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL2_1(:,ii,jj)),ii,0);
        pval2_3_tmp(jj,ii) = test_CR(-log(PL0(:,ii,jj)),-log(PL2_3(:,ii,jj)),ii,0);
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
DM_out = [DM_out;series(series_to_eval)',repmat({'BCVAR'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_1_tmp)];
DM_out = [DM_out;cell(1,size(DM_out,2))];
DM_out = [DM_out;series(series_to_eval)',repmat({'BCTRVAR (opt 3 cov)'},numel(series_to_eval),1),repmat({'AR(1)'},numel(series_to_eval),1),num2cell(pval2_3_tmp)];

xlswrite(f_id,PL_out,'PL');
xlswrite(f_id,DM_out,'DM tests (PL)');
    