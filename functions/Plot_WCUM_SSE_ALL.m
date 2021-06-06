function Plot_WCUM_SSE_ALL(VAR_size,h_list)

addpath('functions')
addpath('data')

this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end

%% Prelims
max_h = 12;
series_to_eval  = 1:7;

%% Prepare data
[Y,series,dates]=Prepare_data(VAR_size);
[T,M] = size(Y);
T_thres = round(0.5*T);

% Covert datees from strings to Matlab dates
dates = datenum(char(dates),'dd/mm/yyyy');

%% Load forecasts from BCVAR and competitors
load([pwd,'/Output/BBENCH_',VAR_size,'.mat']);
load([pwd,'/Output/BVARMINN_',VAR_size,'.mat']);

%No compression in Cov
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_10_1_3_2_0.mat']);
fore2_1     = fore2;
msfe2_1     = msfe2;
msfe2_1_ALL = msfe2_ALL;
PL2_1       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');

%With compression in Cov
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_10_1_3_2_1.mat']);
fore2_3     = fore2;
msfe2_3     = msfe2;
msfe2_3_ALL = msfe2_ALL;
PL2_3       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');

% Compute weighted MSFE

for h=1:size(msfe0,2)
    W = diag(1./diag(cov(Y(T_thres+h:T-max_h+h,:))));
    for t=1:size(msfe0_ALL,1)
        wmsfe0(t,h)   = squeeze(sqrt(msfe0_ALL(t,h,:)))'*W*squeeze(sqrt(msfe0_ALL(t,h,:)));
        wmsfe3(t,h)   = squeeze(sqrt(msfe3_ALL(t,h,:)))'*W*squeeze(sqrt(msfe3_ALL(t,h,:)));
        wmsfe2_1(t,h) = squeeze(sqrt(msfe2_1_ALL(t,h,:)))'*W*squeeze(sqrt(msfe2_1_ALL(t,h,:)));
        wmsfe2_3(t,h) = squeeze(sqrt(msfe2_3_ALL(t,h,:)))'*W*squeeze(sqrt(msfe2_3_ALL(t,h,:)));
        
    end
end

% New figure is needed
fullscreen = get(0,'ScreenSize');
hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);

% WCUM SSE plots
for i=1:length(h_list)
    this_h=h_list(i);
    subplot(ceil(numel(h_list)/2),2,i)
    plot(dates(T_thres+this_h:T-h+this_h),cumsum(-wmsfe3(:,this_h) + wmsfe0(:,this_h)),'linewidth',2)
    hold on;
    plot(dates(T_thres+this_h:T-h+this_h),cumsum(-wmsfe2_1(:,this_h) + wmsfe0(:,this_h)),'linewidth',3)
    plot(dates(T_thres+this_h:T-h+this_h),cumsum(-wmsfe2_3(:,this_h) + wmsfe0(:,this_h)),'linewidth',3)
    plot(dates(T_thres+this_h:T),zeros(T-T_thres-this_h+1,1),'k--','linewidth',2);
    if this_h==1
        legend('BVAR','BCVAR','BCVAR_c','Location','NorthWest');
        legend('boxoff');
    end
    title(['h=',num2str(this_h)])
    ylabel('Weighted Cum SSE');
    xlim([datenum('1-1-1987');datenum('12/1/2014')])
    set(gca,'xtick',[datenum('1-1-1987'):1840:datenum('12/1/2014')]);
    set(gca,'xtick',[datenum('1-1-1987'):1840:datenum('12/1/2014')]);
    datetick('x','yyyy','keepticks','keeplimits')
    set(gca,'Xgrid','on','YGrid','on')
    set(gca,'FontSize',10);
    h_ylabel = get(gca,'YLabel');
    set(h_ylabel,'FontSize',14);
    h_title = get(gca,'Title');
    set(h_title,'FontSize',14);
    
end

% Save figure
f_id_tmp = [this_out,'Figure_',VAR_size,'_WCUM_SSE_ALL'];
set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
%[result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
%delete([f_id_tmp,'.eps']);
close all;
pause(1);



