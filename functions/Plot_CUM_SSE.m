function Plot_CUM_SSE(VAR_size,h_list,date_start,date_end)

addpath('functions')
addpath('data')

this_out  = [pwd,'\Output\',datestr(now,'yyyy.mm.dd'),' Tables and charts\'];
if ~exist(this_out,'dir')
    mkdir(this_out)
end

file_id_all = [VAR_size,'_CUM_SSE_figures.pdf'];
if exist([this_out,'\',file_id_all],'file')
    delete([this_out,'\',file_id_all]);
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

T_beg   = find(dates==datenum(date_start));  
T_end   = find(dates==datenum(date_end));

this_start = T_beg - T_thres + 1;
this_end = T_end - h - T_beg + 1;

%% Load forecasts from BCVAR and competitors
% BENCH
load([pwd,'/Output/BBENCH_',VAR_size,'.mat']);
% DFM and FAVAR models
load([pwd,'/Output/BCOMP_',VAR_size,'.mat']);
% BCVAR
msfe1     = msfe0;
msfe1_ALL = msfe0_ALL;
PL1       = PL0;
fore1     = fore0;

% if ~isequal(VAR_size,'JURADO')
%     load([pwd,'/Output/BCVAR_',VAR_size,'_1_1.mat']);
% else
%     msfe1     = msfe0;
%     msfe1_ALL = msfe0_ALL;
%     PL1       = PL0;
%     fore1     = fore0;
% end

%BCTRVAR, no cov compression, constant coeffs
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_50_1_3_2_0.mat']);
fore2_0     = fore2;
msfe2_0     = msfe2;
msfe2_0_ALL = msfe2_ALL;
PL2_0       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');

%BCTRVAR, cov compression, constant coeffs
load([pwd,'/Output/BCTRVAR_',VAR_size,'_1_50_1_3_2_1.mat']);
fore2_1     = fore2;
msfe2_1     = msfe2;
msfe2_1_ALL = msfe2_ALL;
PL2_1       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');

%BCTRVAR, with cov compression and TVP-SV
load([pwd,'/Output/BCTRVAR_TVP_',VAR_size,'_1_50_1_3_2_1.mat']);
fore2_3     = fore2;
msfe2_3     = msfe2;
msfe2_3_ALL = msfe2_ALL;
PL2_3       = PL2;
clear('fore2','msfe2','msfe2_ALL','PL2');
% BVAR with MINNESOTA priors
load([pwd,'/Output/BVARMINN_',VAR_size,'.mat']);


% CUM SSE plots
for this_h=h_list
    % New figure is needed
    fullscreen = get(0,'ScreenSize');
    hh = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
    for i=1:length(series_to_eval)
    
        subplot(ceil(numel(series_to_eval)/2),2,i)
        plot(dates(T_beg+this_h:T_end-h+this_h),cumsum(-msfe4(this_start:this_end,this_h,i) + msfe0(this_start:this_end,this_h,i)),'linewidth',2)
        hold on;
        plot(dates(T_beg+this_h:T_end-h+this_h),cumsum(-msfe5(this_start:this_end,this_h,i) + msfe0(this_start:this_end,this_h,i)),'linewidth',2)
        plot(dates(T_beg+this_h:T_end-h+this_h),cumsum(-msfe3(this_start:this_end,this_h,i) + msfe0(this_start:this_end,this_h,i)),'linewidth',2)
        plot(dates(T_beg+this_h:T_end-h+this_h),cumsum(-msfe2_0(this_start:this_end,this_h,i) + msfe0(this_start:this_end,this_h,i)),'linewidth',2)
        plot(dates(T_beg+this_h:T_end-h+this_h),cumsum(-msfe2_1(this_start:this_end,this_h,i) + msfe0(this_start:this_end,this_h,i)),'linewidth',2)
        plot(dates(T_beg+this_h:T_end-h+this_h),cumsum(-msfe2_3(this_start:this_end,this_h,i) + msfe0(this_start:this_end,this_h,i)),'linewidth',2)
        
        plot(dates(T_beg+this_h:T_end-h+this_h),zeros(T_end-T_beg-h+1,1),'k--','linewidth',2);
        if i==1
            legend('DFM','FAVAR','BVAR','BCVAR','BCVAR_c','BCVAR_{tvp}','Location','NorthWest','Orientation','horizontal');
            legend('boxoff');
        end
        title([char(series(series_to_eval(i))),' h=',num2str(this_h)])
        ylabel('CSSFED');
        xlim([datenum(dates(T_beg+this_h));datenum(dates(T_end-h+this_h))])
        set(gca,'xtick',[datenum(dates(T_beg+this_h)):1840:datenum(dates(T_end-h+this_h))]);
        set(gca,'xtick',[datenum(dates(T_beg+this_h)):1840:datenum(dates(T_end-h+this_h))]);
        datetick('x','yyyy','keepticks','keeplimits')
        set(gca,'Xgrid','on','YGrid','on')
        set(gca,'FontSize',10);
        h_ylabel = get(gca,'YLabel');
        set(h_ylabel,'FontSize',14);
        h_title = get(gca,'Title');
        set(h_title,'FontSize',14);
    end
    
    % Save figure
    f_id_tmp = [this_out,'Figure ',VAR_size,' CUM SSE h=',num2str(this_h),' (',datestr(date_start,'yyyy'),'--',datestr(date_end,'yyyy'),')'];
    set(gcf,'PaperPositionMode','auto')
    print('-depsc','-tiff','-r600',[f_id_tmp,'.eps'])
    [result,msg] = eps2xxx([f_id_tmp,'.eps'],{'pdf'},'C:\Program Files\gs\gs9.04\bin\gswin64c.exe');
    delete([f_id_tmp,'.eps']);
    close all;
    pause(1);
    append_pdfs([this_out,'\',file_id_all], [f_id_tmp,'.pdf'])
    delete([f_id_tmp,'.pdf']);
end


    