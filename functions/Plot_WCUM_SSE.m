function Plot_WCUM_SSE(VAR_size,h_list,date_start,date_end)

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

T_beg   = find(dates==datenum(date_start));  
T_end   = find(dates==datenum(date_end));

this_start = T_beg - T_thres + 1;
this_end = T_end - max_h - T_beg + 1;


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





% Compute weighted MSFE

for h=1:size(msfe0,2)
    W = diag(1./diag(cov(Y(T_beg+h:T_end-max_h+h,series_to_eval))));
    for t=this_start:this_end
        wmsfe0(t,h)   = squeeze(sqrt(msfe0(t,h,:)))'*W*squeeze(sqrt(msfe0(t,h,:)));
        wmsfe1(t,h)   = squeeze(sqrt(msfe1(t,h,:)))'*W*squeeze(sqrt(msfe1(t,h,:)));
        wmsfe3(t,h)   = squeeze(sqrt(msfe3(t,h,:)))'*W*squeeze(sqrt(msfe3(t,h,:)));
        wmsfe4(t,h)   = squeeze(sqrt(msfe4(t,h,:)))'*W*squeeze(sqrt(msfe4(t,h,:)));
        wmsfe5(t,h)   = squeeze(sqrt(msfe5(t,h,:)))'*W*squeeze(sqrt(msfe5(t,h,:)));
        wmsfe6(t,h)   = squeeze(sqrt(msfe6(t,h,:)))'*W*squeeze(sqrt(msfe6(t,h,:)));
        wmsfe2_1(t,h) = squeeze(sqrt(msfe2_1(t,h,:)))'*W*squeeze(sqrt(msfe2_1(t,h,:)));
        wmsfe2_3(t,h) = squeeze(sqrt(msfe2_3(t,h,:)))'*W*squeeze(sqrt(msfe2_3(t,h,:)));
        
    end
end

% New figure is needed
fullscreen = get(0,'ScreenSize');
hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);

% WCUM SSE plots
for i=1:length(h_list)
    this_h=h_list(i);
    subplot(ceil(numel(h_list)/2),2,i)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(-wmsfe4(this_start:this_end,this_h) + wmsfe0(this_start:this_end,this_h)),'linewidth',2)
    hold on;
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(-wmsfe5(this_start:this_end,this_h) + wmsfe0(this_start:this_end,this_h)),'linewidth',2)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(-wmsfe3(this_start:this_end,this_h) + wmsfe0(this_start:this_end,this_h)),'linewidth',2)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(-wmsfe1(this_start:this_end,this_h) + wmsfe0(this_start:this_end,this_h)),'linewidth',3)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(-wmsfe2_1(this_start:this_end,this_h) + wmsfe0(this_start:this_end,this_h)),'linewidth',3)
    %plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(-wmsfe2_3(this_start:this_end,this_h) + wmsfe0(this_start:this_end,this_h)),'linewidth',3)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),zeros(T_end-T_beg-max_h+1,1),'k--','linewidth',2);
    if this_h==1
        %legend('DFM','FAVAR','BVAR','BCVAR','BCVAR_c','BCVAR_{tvp}','Location','SouthWest','Orientation','horizontal');
        legend('DFM','FAVAR','BVAR','BCVAR','BCVAR_c','Location','SouthWest','Orientation','horizontal');
        legend('boxoff');
    end
    title(['h=',num2str(this_h)])
    ylabel('CSWFED');
    
    xlim([datenum(dates(T_beg+this_h));datenum(dates(T_end-max_h+this_h))])
    set(gca,'xtick',[datenum(dates(T_beg+this_h)):1840:datenum(dates(T_end-max_h+this_h))]);
    set(gca,'xtick',[datenum(dates(T_beg+this_h)):1840:datenum(dates(T_end-max_h+this_h))]);
    datetick('x','yyyy','keepticks','keeplimits')
    set(gca,'Xgrid','on','YGrid','on')
    set(gca,'FontSize',10);
    h_ylabel = get(gca,'YLabel');
    set(h_ylabel,'FontSize',14);
    h_title = get(gca,'Title');
    set(h_title,'FontSize',14);
    
end

% Save figure
f_id_tmp = [this_out,'Figure_',VAR_size,'_WCUM_SSE_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy')];
set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
[result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
delete([f_id_tmp,'.eps']);
close all;
pause(1);



