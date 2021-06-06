% PURPOSE: Evaluate PL from competing models
clear;
clc;

addpath('functions')
addpath('data')

%% Preliminaries
VAR_size = 'MEDIUM';
h = 12;
p = 13;
ndraws = 1000;
series_to_eval = 1:7;

% settings for BCTRVAR
RP_type       = 1;
n_psi         = 50; 
stdata        = 1;  % 0: do nothing; 1: standardize data;

apply_bcr     = 3;  % 1:everywhere; 2:intercepts excluded; 3: intercepts and first own lags excluded
weight_scheme = 2;  % 1: 1/N; 2: BIC for whole VAR; 3: equation-by-equation BIC
cov_comp      = 0;  % 0: Don't include covariance matrix terms; 1: Include cov matrix terms
sparsity      = 0;

this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\Kdensities'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end

%% Prepare data
[Y,series,dates]=Prepare_data(VAR_size);
[T,M] = size(Y);
T_thres = round(0.5*T);


% Load output mats, PL0 is benchmark, PL2 is from BCTRVAR
load([pwd,'\Output\BBENCH_VAR_',VAR_size,'.mat']);
% load([pwd,'\Output\BCTRVAR_',VAR_size,'_1_10_1_3_2_0 - 1 to 5ln(N).mat']);
% PL2_0   = PL2;
% fore2_0 = fore2;
% clear('PL2','fore2');
load([pwd,'\Output\BCTRVAR_',VAR_size,'_1_10_1_3_2_0.mat']);
PL2_1   = PL2;
fore2_1 = fore2;
clear('PL2','fore2');
load([pwd,'\Output\BVARMINN_',VAR_size,'.mat']);
load([pwd,'\Output\BCOMP_',VAR_size,'.mat']);

% Some output

for this_h=1  %1:h
    dates_out  = datenum(dates(T_thres+this_h:end-h+this_h),'dd/mm/yyyy');
    
    this_PL0   = squeeze(PL0(:,this_h,:));
    %this_PL2_0 = squeeze(PL2_0(:,this_h,:));
    this_PL2_1 = squeeze(PL2_1(:,this_h,:));
    this_PL3   = squeeze(PL3(:,this_h,:));
    this_PL5   = squeeze(PL5(:,this_h,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    this_PL0(this_PL0<1e-2)     = 1e-2;
    %this_PL2_0(this_PL2_0<1e-2) = 1e-2;
    this_PL2_1(this_PL2_1<1e-2) = 1e-2;
    this_PL3(this_PL3<1e-2)     = 1e-2;
    this_PL5(this_PL5<1e-2)     = 1e-2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculating log predictive bayes factors; note: these are
    % relative to the corresponding benchmark model;
    %log_BF2_0   = log(this_PL2_0./this_PL0);
    log_BF2_1   = log(this_PL2_1./this_PL0);
    log_BF3     = log(this_PL3./this_PL0);
    log_BF5     = log(this_PL5./this_PL0);
    
    % Set to zero all cases when Predictive likelihood function is zero
    % for no predictability model (to avoid Inf from the ratios)
    %log_BF2_0(this_PL0==0) = 0;
    log_BF2_1(this_PL0==0) = 0;
    log_BF3(this_PL0==0)   = 0;
    log_BF5(this_PL0==0)   = 0;
    
    % Set NaNs to zero, just for plotting purposes
    %log_BF2_0(isnan(log_BF2_0)) = 0;
    log_BF2_1(isnan(log_BF2_1)) = 0;
    log_BF3(isnan(log_BF3)) = 0;
    log_BF5(isnan(log_BF3)) = 0;
    
    for t=T_thres:T-h
        
        yi = Y(t+this_h,series_to_eval);
        
        this_fore0   = squeeze(fore0(t-T_thres+1,:,this_h,:));
        %this_fore2_0 = squeeze(fore2_0(t-T_thres+1,:,this_h,:));
        this_fore2_1 = squeeze(fore2_1(t-T_thres+1,:,this_h,:));
        this_fore3   = squeeze(fore3(t-T_thres+1,:,this_h,:));
        this_fore5   = squeeze(fore5(t-T_thres+1,:,this_h,:));
    
        % New figure is needed
        fullscreen = get(0,'ScreenSize');
        hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
        
        for j=1:length(series_to_eval)
            this_series = series(j);
            subplot(ceil(length(series_to_eval)/2),2,j);
            
            ksdensity(this_fore0(:,j));
            hold on
            %ksdensity(this_fore2_0(:,j));
            ksdensity(this_fore2_1(:,j));
            ksdensity(this_fore3(:,j));
            ksdensity(this_fore5(:,j));
            max_y = ylim;
            stem(yi(j),max_y(2),'filled')
            title([datestr(dates_out(t-T_thres+1),'mmm-yyyy'),' -- ',char(this_series),' h=',num2str(this_h)])
            ylabel('Rel. freq.');
            if j==1
                legend('VAR OLS','BCVAR','BVAR','FAVAR','location','NorthWest');
                %legend('boxoff');
            end
            set(gca,'Xgrid','on','YGrid','on')
            set(gca,'FontSize',10);
            h_ylabel = get(gca,'YLabel');
            set(h_ylabel,'FontSize',14);
            h_title = get(gca,'Title');
            set(h_title,'FontSize',14);
        end
    
        % Save figure
        f_id_tmp = [this_out,'\',datestr(dates_out(t-T_thres+1),'yyyy.mm'),' Kdensity ',VAR_size,' h=',num2str(this_h)];
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
        [result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
        delete([f_id_tmp,'.eps']);
        close all;

    end
end