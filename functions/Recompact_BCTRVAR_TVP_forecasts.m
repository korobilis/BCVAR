function Recompact_BCTRVAR_TVP_forecasts(VAR_size,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp)

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
    f_id = [pwd,'/Temp/',sprintf('%s_%s_%g_%g_%g_%g_%g_%g_%s.mat','BCTRVAR_TVP',VAR_size,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp,irep_to_print)];
    load(f_id);
    
    for ii = 1:h
        % Save forecasts
        fore2(irep-T_thres+1,:,ii,:) = bcvarcon2(:,ii,series_to_eval);
        
        yi = Y(irep+ii,series_to_eval);
        
        msfe2(irep-T_thres+1,ii,:) = (squeeze(mean(bcvarcon2(:,ii,series_to_eval),1))' - yi).^2;
        msfe2_ALL(irep-T_thres+1,ii,:) = (squeeze(mean(bcvarcon2(:,ii,:),1))' - Y(irep+ii,:)).^2;
        
        if size(bcvarcon2,1) > 1
            for j = series_to_eval
                PL2(irep-T_thres+1,ii,j) = ksdensity(squeeze(fore2(irep-T_thres+1,:,ii,j)),yi(:,j));
            end
        else
            PL2(irep-T_thres+1,ii,series_to_eval) = NaN;
        end
    end
    clear('bcvarcon2');
    %delete(f_id);
    
end
    
%% Save output
save([pwd,'/Output/',sprintf('%s_%s_%g_%g_%g_%g_%g_%g.mat','BCTRVAR_TVP',VAR_size,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp)],'Y','fore*','msfe*','PL*','-mat');
