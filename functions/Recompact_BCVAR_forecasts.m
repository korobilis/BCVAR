function Recompact_BCVAR_forecasts(VAR_size,RP_type,n_psi)

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
    f_id = [pwd,'/Temp/',sprintf('%s_%s_%g_%g_%s.mat','BCVAR',VAR_size,RP_type,n_psi,irep_to_print)];
    load(f_id);
    
    for ii = 1:h
        % Save forecasts
        fore1(irep-T_thres+1,:,ii,:) = bcvarcon1(:,ii,series_to_eval);
        
        yi = Y(irep+ii,series_to_eval);
        
        msfe1(irep-T_thres+1,ii,:) = (squeeze(mean(bcvarcon1(:,ii,series_to_eval),1))' - yi).^2;
        msfe1_ALL(irep-T_thres+1,ii,:) = (squeeze(mean(bcvarcon1(:,ii,:),1))' - Y(irep+ii,:)).^2;
        
        if size(bcvarcon1,1) > 1
            for j = series_to_eval
                PL1(irep-T_thres+1,ii,j) = ksdensity(squeeze(fore1(irep-T_thres+1,:,ii,j)),yi(:,j));
            end
        else
            PL1(irep-T_thres+1,ii,series_to_eval) = NaN;
        end
    end
    clear('bcvarcon1');
    delete(f_id);
    
end
    
%% Save output
save([pwd,'/Output/',sprintf('%s_%s_%g_%g.mat','BCVAR',VAR_size,RP_type,n_psi)],'Y','fore*','msfe*','PL*','-mat');
