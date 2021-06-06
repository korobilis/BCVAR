function Plot_forecasts(VAR_size,h_list)

addpath('functions')
addpath('data')

this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end

%% Prelims
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
load([pwd,'/Output/BCTRVAR_TVP_',VAR_size,'_1_10_1_3_2_1.mat']);
fore2_3     = fore2;
msfe2_3     = msfe2;
msfe2_3_ALL = msfe2_ALL;
PL2_3       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');
% BVAR with MINNESOTA priors
load([pwd,'/Output/BVARMINN_',VAR_size,'.mat']);

% Prep LHS series
for irep = T_thres:T-h
    for ii = 1:h
        yi(irep-T_thres+1,:,ii) = Y(irep+ii,series_to_eval);
    end
end
keyboard
% Forecast plots
for this_h=h_list
    %AR(1) benchmark
    if ndims(fore0) == 4
        this_forc(:,:,1) = squeeze(mean(squeeze(fore0(:,:,this_h,:)),2));
        this_std(:,:,1)  = squeeze(std(squeeze(fore0(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,1) = fore0(:,:,this_h,:);
    end
    % DFM
    if ndims(fore4) == 4
        this_forc(:,:,2) = squeeze(mean(squeeze(fore4(:,:,this_h,:)),2));
        this_std(:,:,2)  = squeeze(std(squeeze(fore4(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,2) = fore4(:,:,this_h,:);
    end
    % FAVAR1
    if ndims(fore5) == 4
        this_forc(:,:,3) = squeeze(mean(squeeze(fore5(:,:,this_h,:)),2));
        this_std(:,:,3)  = squeeze(std(squeeze(fore5(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,3) = fore5(:,:,this_h,:);
    end
    %BVAR-MINN
    if ndims(fore3) == 4
        this_forc(:,:,4) = squeeze(mean(squeeze(fore3(:,:,this_h,:)),2));
        this_std(:,:,4)  = squeeze(std(squeeze(fore3(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,4) = fore3(:,:,this_h);
    end
    %BCTRVAR (option 3, no cov)
    if ndims(fore1) == 4
        this_forc(:,:,5) = squeeze(mean(squeeze(fore1(:,:,this_h,:)),2));
        this_std(:,:,5)  = squeeze(std(squeeze(fore1(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,5)        = fore1(:,:,this_h,:);
    end
    %BCTRVAR (option 3, cov)
    if ndims(fore2_1) == 4
        this_forc(:,:,6) = squeeze(mean(squeeze(fore2_1(:,:,this_h,:)),2));
        this_std(:,:,6)  = squeeze(std(squeeze(fore2_1(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,6)        = fore2_1(:,:,this_h,:);
    end
    %BCTRVAR TVP (option 3, cov)
    if ndims(fore2_3) == 4
        this_forc(:,:,7) = squeeze(mean(squeeze(fore2_3(:,:,this_h,:)),2));
        this_std(:,:,7)  = squeeze(std(squeeze(fore2_3(:,:,this_h,:)),0,2));
    else
        this_forc(:,:,7)        = fore2_3(:,:,this_h,:);
    end
%     % New figure is needed
%     fullscreen = get(0,'ScreenSize');
%     hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
%     for i=1:length(series_to_eval)
%         subplot(ceil(numel(series_to_eval)/2),2,i)
%         plot(dates(T_thres+this_h:T-h+this_h),yi(:,i,this_h),'.k','MarkerSize',20);
%         hold on;
%         h1 = plot(dates(T_thres+this_h:T-h+this_h),yi(:,i,this_h),'k');
%         h2 = plot(dates(T_thres+this_h:T-h+this_h),this_forc(:,i,3),'linewidth',2);
%         h3 = plot(dates(T_thres+this_h:T-h+this_h),this_forc(:,i,4),'linewidth',2);
%         h4 = plot(dates(T_thres+this_h:T-h+this_h),this_forc(:,i,6),'linewidth',2);
%         h5 = plot(dates(T_thres+this_h:T-h+this_h),this_forc(:,i,7),'linewidth',2);
%         
%         if i==1
%             legend([h1,h2,h3,h4,h5],'Actual','FAVAR1','BVAR MINN','BCVAR_c','BCVAR_{tvp}','Location','SouthWest','Orientation','Horizontal');
%             legend('boxoff');
%         end
%         title([char(series(series_to_eval(i))),' - h=',num2str(this_h)])
%         xlim([datenum('1-1-1982');datenum('12/1/2014')])
%         set(gca,'xtick',[datenum('1-1-1982'):1840:datenum('12/1/2014')]);
%         set(gca,'xtick',[datenum('1-1-1982'):1840:datenum('12/1/2014')]);
%         datetick('x','yyyy','keepticks','keeplimits')
%         set(gca,'Xgrid','on','YGrid','on')
%         set(gca,'FontSize',10);
%         h_ylabel = get(gca,'YLabel');
%         set(h_ylabel,'FontSize',14);
%         h_title = get(gca,'Title');
%         set(h_title,'FontSize',14);
%        
%     end
%     % Save figure
%     f_id_tmp = [this_out,'Figure ',VAR_size,' forc h=',num2str(this_h)];
%     set(gcf,'PaperPositionMode','auto')
%     print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
%     [result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
%     delete([f_id_tmp,'.eps']);
%     close all;
%     pause(1);
    
    % New figure is needed
    fullscreen = get(0,'ScreenSize');
    hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
    for i=1:length(series_to_eval)
        subplot(ceil(numel(series_to_eval)/2),2,i)
        h1 = plot(dates(T_thres+this_h:T-h+this_h),this_std(:,i,2),'linewidth',2);
        hold on;
        h2 = plot(dates(T_thres+this_h:T-h+this_h),this_std(:,i,3),'linewidth',2);
        h3 = plot(dates(T_thres+this_h:T-h+this_h),this_std(:,i,4),'linewidth',2);
        h4 = plot(dates(T_thres+this_h:T-h+this_h),this_std(:,i,5),'linewidth',2);
        h5 = plot(dates(T_thres+this_h:T-h+this_h),this_std(:,i,6),'linewidth',2);
        h6 = plot(dates(T_thres+this_h:T-h+this_h),this_std(:,i,7),'linewidth',2);
        
        if i==1
            legend([h1,h2,h3,h4,h5,h6],'DFM','FAVAR','BVAR','BCVAR','BCVAR_c','BCVAR_{tvp}','Location','NorthWest','Orientation','Horizontal');
            legend('boxoff');
        end
        title([char(series(series_to_eval(i))),' - h=',num2str(this_h)])
        xlim([datenum('1-1-1982');datenum('12/1/2014')])
        set(gca,'xtick',[datenum('1-1-1982'):1840:datenum('12/1/2014')]);
        set(gca,'xtick',[datenum('1-1-1982'):1840:datenum('12/1/2014')]);
        datetick('x','yyyy','keepticks','keeplimits')
        set(gca,'Xgrid','on','YGrid','on')
        set(gca,'FontSize',10);
        h_ylabel = get(gca,'YLabel');
        set(h_ylabel,'FontSize',14);
        h_title = get(gca,'Title');
        set(h_title,'FontSize',14);
       
    end
end
    