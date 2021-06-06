function Recompact_BVAR_MINN_forecasts(VAR_size)

addpath('functions')
addpath('data')

if ~exist([pwd,'/Output'],'dir')
    mkdir([pwd,'/Output'])
end

[Y,series,dates]=Prepare_data(VAR_size);
[T,M] = size(Y);
T_thres = round(0.5*T);

%% Load the individual forecast files and combine them

%%%%%%
h=12;
%%%%%%

series_to_eval = 1:7;

for irep = T_thres:T-h
    jj = irep-T_thres+1;
    if jj < 10
        irep_to_print = ['00',num2str(jj)];
    elseif jj < 100
        irep_to_print = ['0',num2str(jj)];
    else
        irep_to_print = num2str(jj);
    end
    disp(['Iteration ',num2str(jj),' of ',num2str(T-h-T_thres+1)]);
    f_id = [pwd,'/Temp/',sprintf('%s_%s_%s.mat','BVAR_MINN',VAR_size,irep_to_print)];
    load(f_id);
    
    for ii = 1:h
        % Save forecasts
        fore3(irep-T_thres+1,:,ii,:) = bvarmin(:,ii,series_to_eval);
        
        yi = Y(irep+ii,series_to_eval);
        
        msfe3(irep-T_thres+1,ii,:) = (squeeze(mean(bvarmin(:,ii,series_to_eval),1))' - yi).^2;
        msfe3_ALL(irep-T_thres+1,ii,:) = (squeeze(mean(bvarmin(:,ii,:),1))' - Y(irep+ii,:)).^2;
        
        if size(bvarmin,1) > 1
            for j = series_to_eval
                PL3(irep-T_thres+1,ii,j) = ksdensity(squeeze(fore3(irep-T_thres+1,:,ii,j)),yi(:,j));
            end
        else
            PL3(irep-T_thres+1,ii,series_to_eval) = NaN;
        end
            
    end
    clear('bvarmin');
    %delete(f_id);
    
end
    
%% Save output
save([pwd,'/Output/',sprintf('%s_%s.mat','BVARMINN',VAR_size)],'Y','fore*','msfe*','PL*','-mat');
