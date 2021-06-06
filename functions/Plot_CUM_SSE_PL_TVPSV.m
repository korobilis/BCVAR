function Plot_CUM_SSE_PL_TVPSV(VAR_size_list,date_start,date_end)

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
[Y,series,dates]=Prepare_data(char(VAR_size_list(1)));
[T,M] = size(Y);
T_thres = round(0.5*T);

% Covert datees from strings to Matlab dates
dates = datenum(char(dates),'dd/mm/yyyy');

T_beg   = find(dates==datenum(date_start));  
T_end   = find(dates==datenum(date_end));

this_start = T_beg - T_thres + 1;
this_end = T_end - max_h - T_beg + 1;
% New figure is needed
fullscreen = get(0,'ScreenSize');
hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);

for i=1:length(VAR_size_list)
    VAR_size = char(VAR_size_list(i));

    %% Load forecasts from BCVAR and competitors
    % BENCH
    load([pwd,'/Output/BBENCH_',VAR_size,'.mat']);

    %BCTRVAR with TVP and SV, with cov compression
    load([pwd,'/Output/BCTRVAR_TVP_',VAR_size,'_1_10_1_3_2_1.mat']);
    fore2_3     = fore2;
    msfe2_3     = msfe2;
    msfe2_3_ALL = msfe2_ALL;
    PL2_3       = PL2;
    clear('fore2','msfe2','msfe2_ALL','PL2');
    
    % Compute weighted MSFE

    for h=1:size(msfe0,2)
        W = diag(1./diag(cov(Y(T_beg+h:T_end-max_h+h,series_to_eval))));
        for t=this_start:this_end
            wmsfe0(t,h)   = squeeze(sqrt(msfe0(t,h,:)))'*W*squeeze(sqrt(msfe0(t,h,:)));
            wmsfe2_3(t,h) = squeeze(sqrt(msfe2_3(t,h,:)))'*W*squeeze(sqrt(msfe2_3(t,h,:)));

        end
    end
    
    % Load Multivariate PLs
    load([pwd,'/Output/MV-PLS_',VAR_size]);
    

    % WCUM SSE and CUM WPL plots
    subplot(numel(VAR_size_list),2,(i-1)*2+1)
    plot(dates(T_beg+1:T_end-max_h+1),cumsum(-wmsfe2_3(this_start:this_end,1) + wmsfe0(this_start:this_end,1)),'linewidth',3)
    hold on;
    plot(dates(T_beg+12:T_end-max_h+12),cumsum(-wmsfe2_3(this_start:this_end,12) + wmsfe0(this_start:this_end,12)),'linewidth',3)
    plot(dates(T_beg:T_end),zeros(T_end-T_beg+1,1),'k--','linewidth',2);
    if i==1
        legend('h=1','h=12','Location','NorthWest','Orientation','vertical');
        legend('boxoff');
    end
    title(['Weighted Cum SSE differentials - ',VAR_size,' VAR'])
    ylabel('Weighted Cum SSE');
    xlim([datenum(dates(T_beg+1));datenum(dates(T_end-max_h+12))])
    set(gca,'xtick',[datenum(dates(T_beg+1)):1840:datenum(dates(T_end-max_h+12))]);
    set(gca,'xtick',[datenum(dates(T_beg+1)):1840:datenum(dates(T_end-max_h+12))]);
    datetick('x','yyyy','keepticks','keeplimits')
    set(gca,'Xgrid','on','YGrid','on')
    set(gca,'FontSize',10);
    h_ylabel = get(gca,'YLabel');
    set(h_ylabel,'FontSize',14);
    h_title = get(gca,'Title');
    set(h_title,'FontSize',14);
    
    subplot(numel(VAR_size_list),2,i*2)
    plot(dates(T_beg+1:T_end-max_h+1),cumsum(log(PL2_3_all(this_start:this_end,1)) - log(PL0_all(this_start:this_end,1))),'linewidth',3)
    hold on;
    plot(dates(T_beg+12:T_end-max_h+12),cumsum(log(PL2_3_all(this_start:this_end,12)) - log(PL0_all(this_start:this_end,12))),'linewidth',3)
    plot(dates(T_beg:T_end),zeros(T_end-T_beg+1,1),'k--','linewidth',2);
    title(['Multivariate Cum PL differentials - ',VAR_size,' VAR'])
    ylabel('Multivariate Cum PL');
    xlim([datenum(dates(T_beg+1));datenum(dates(T_end-max_h+12))])
    set(gca,'xtick',[datenum(dates(T_beg+1)):1840:datenum(dates(T_end-max_h+12))]);
    set(gca,'xtick',[datenum(dates(T_beg+1)):1840:datenum(dates(T_end-max_h+12))]);
    datetick('x','yyyy','keepticks','keeplimits')
    set(gca,'Xgrid','on','YGrid','on')
    set(gca,'FontSize',10);
    h_ylabel = get(gca,'YLabel');
    set(h_ylabel,'FontSize',14);
    h_title = get(gca,'Title');
    set(h_title,'FontSize',14);
    
end

% Save figure
f_id_tmp = [this_out,'Figure_BCVAR_TVPSV_',datestr(date_start,'yyyy'),'_',datestr(date_end,'yyyy')];
set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
[result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
delete([f_id_tmp,'.eps']);
close all;
pause(1);



