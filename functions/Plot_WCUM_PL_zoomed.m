function Plot_WCUM_PL_zoomed(VAR_size,h_list,date_start,date_end)

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

%% Load Multivariate PLs
load([pwd,'/Output/MV-PLS_',VAR_size]);

% New figure is needed
fullscreen = get(0,'ScreenSize');
hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);

% WCUM SSE plots
for i=1:length(h_list)
    this_h=h_list(i);
    subplot(ceil(numel(h_list)/2),2,i)
    
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(log(PL4_all(this_start:this_end,this_h)) - log(PL0_all(this_start:this_end,this_h))),'linewidth',2)
    hold on;
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(log(PL5_all(this_start:this_end,this_h)) - log(PL0_all(this_start:this_end,this_h))),'linewidth',2)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(log(PL3_all(this_start:this_end,this_h)) - log(PL0_all(this_start:this_end,this_h))),'linewidth',2)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(log(PL1_all(this_start:this_end,this_h)) - log(PL0_all(this_start:this_end,this_h))),'linewidth',3)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(log(PL2_1_all(this_start:this_end,this_h)) - log(PL0_all(this_start:this_end,this_h))),'linewidth',3)
    %plot(dates(T_beg+this_h:T_end-max_h+this_h),cumsum(log(PL2_3_all(this_start:this_end,this_h)) - log(PL0_all(this_start:this_end,this_h))),'linewidth',3)
    plot(dates(T_beg+this_h:T_end-max_h+this_h),zeros(T_end-T_beg-max_h+1,1),'k--','linewidth',2);
    if this_h==1
        %legend('DFM','FAVAR','BVAR','BCVAR','BCVAR_c','BCVAR_{tvp}','Location','NorthWest','Orientation','horizontal');
        legend('DFM','FAVAR','BVAR','BCVAR','BCVAR_c','Location','NorthWest','Orientation','horizontal');
        legend('boxoff');
    end
    title(['h=',num2str(this_h)])
    ylabel('CSMVLPLD');
    xlim([datenum('1/1/2007');datenum('1/12/2012')])
    set(gca,'xtick',[datenum('1/1/2007'):90:datenum('1/12/2012')]);
    datetick('x','mmm-yyyy','keepticks','keeplimits')
    rotateXLabels( gca, 45)
    set(gca,'Xgrid','on','YGrid','on')
    set(gca,'FontSize',10);
    h_ylabel = get(gca,'YLabel');
    set(h_ylabel,'FontSize',14);
    h_title = get(gca,'Title');
    set(h_title,'FontSize',14);
    
end

% Save figure
f_id_tmp = [this_out,'Figure_',VAR_size,'_WCUM_PL_zoomed_in_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy')];
set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
[result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
delete([f_id_tmp,'.eps']);
close all;
pause(1);



