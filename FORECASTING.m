% FORECASTING.m    This code demonstrates the recursive forecasts from the Compressed VAR and alternative benchmark
% models using in Koop, Korobilis, Pettenuzzo (2017) Bayesian compressed vector autoregressions, Journal of Econometrics
% ======================================================================================================================

clear;
clc;

addpath('functions')
addpath('data')

%% Preliminaries
VAR_size = 'MEDIUM';
h = 12;
p = 13;
ndraws = 1;    % ***** ATTENTION: Typically you would want around 10000 draws from the predictive density, but don't set 
               % a large value unless you know you have many GBs or available RAM****
               % If you set ndraws = 1, then the code automatically does forecasting using the analytically available
               % posterior predictive mean for each competing method
series_to_eval = 1:7;


%% Prepare data
[Y,series,dates]=Prepare_data(VAR_size);
[T,M] = size(Y);
T_thres = round(0.5*T);

% settings for VAR models
RP_type       = 1;  % Projection type (see genRP.m)
n_psi         = 50; % Number of Random Projections
stdata        = 1;  % 0: do nothing; 1: standardize data;
apply_bcr     = 3;  % 1:everywhere; 2:intercepts excluded; 3: intercepts and first own lags excluded
weight_scheme = 2;  % 1: 1/N; 2: BIC for whole VAR; 3: equation-by-equation BIC
cov_comp      = 0;  % 0: Don't include covariance matrix terms in compression; 1: Include cov matrix terms
sparsity      = 0;

lambda        = 0.99;  % Forgetting factor for TVP-VAR model
delta         = 0.94;  % Decay factor for TVP-VAR model
forgetting    = 1;     % 0: choose lambda,delta; 1: estimate lambda,delta optimally

grid = [0.5:.1:10 50 100];  % grid for Minnesota prior BVAR
lgrid = grid*sqrt(M*p);
maxfac = 2*round(sqrt(M));    % Max number of factors for FAVAR/DFM


fore     = zeros(T-T_thres+1,ndraws,h,length(series_to_eval),8);
msfe     = zeros(T-T_thres+1,h,length(series_to_eval),8);
msfe_ALL = zeros(T-T_thres+1,h,M,8);
PL       = zeros(T-T_thres+1,h,length(series_to_eval),8);

%% Run estimation and forecasts
for irep = T_thres:T-h       
    disp(['Iteration ',num2str(irep-T_thres+1),' of ',num2str(T-h-T_thres+1)]);
    bctrvar     = BCTRVAR_CONJ(Y(1:irep,:),p,1,h,ndraws,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp,sparsity,series_to_eval);
    bctvpvar    = BCTRVAR_TVP(Y(1:irep,:),p,1,h,ndraws,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp,sparsity,series_to_eval,lambda,delta,forgetting);
    bvarmin     = BVAR_MINN(Y(1:irep,:),p,1,h,ndraws,lgrid,0);
    bfavar1     = BFAVAR(Y(1:irep,:),p,1,h,ndraws,series_to_eval,maxfac);
    bfavar2     = BFAVAR(Y(1:irep,:),p,1,h,ndraws,series_to_eval,1);
    bvardfm     = BDFM(Y(1:irep,:),p,1,h,ndraws,maxfac,series_to_eval);
    bvarbnch    = BVAR_OLS(Y(1:irep,series_to_eval),p,h,ndraws); 
    barbnch     = BAR(Y(1:irep,:),p,h,ndraws);
    
    for ii = 1:h
        % Save forecasts
        fore(irep-T_thres+1,:,ii,:,1) = bctrvar(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,2) = bctvpvar(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,3) = bvarmin(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,4) = bfavar1(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,5) = bfavar2(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,6) = bvardfm(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,7) = bvarbnch(:,ii,series_to_eval);
        fore(irep-T_thres+1,:,ii,:,8) = barbnch(:,ii,series_to_eval);
            
        msfe(irep-T_thres+1,ii,:,1)     = (squeeze(mean(bctrvar(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,2)     = (squeeze(mean(bctvpvar(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,3)     = (squeeze(mean(bvarmin(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,4)     = (squeeze(mean(bfavar1(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,5)     = (squeeze(mean(bfavar2(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,6)     = (squeeze(mean(bvardfm(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,7)     = (squeeze(mean(bvarbnch(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        msfe(irep-T_thres+1,ii,:,8)     = (squeeze(mean(barbnch(:,ii,series_to_eval),1))' - Y(irep+ii,series_to_eval)).^2;
        
        msfe_ALL(irep-T_thres+1,ii,:,1) = (squeeze(mean(bctrvar(:,ii,:),1))' - Y(irep+ii,:)).^2;
        msfe_ALL(irep-T_thres+1,ii,:,2) = (squeeze(mean(bctvpvar(:,ii,:),1))' - Y(irep+ii,:)).^2;
        msfe_ALL(irep-T_thres+1,ii,:,3) = (squeeze(mean(bvarmin(:,ii,:),1))' - Y(irep+ii,:)).^2;
%         msfe_ALL(irep-T_thres+1,ii,:,4) = (squeeze(mean(bfavar1(:,ii,:),1))' - Y(irep+ii,:)).^2;
%         msfe_ALL(irep-T_thres+1,ii,:,5) = (squeeze(mean(bfavar2(:,ii,:),1))' - Y(irep+ii,:)).^2;
        msfe_ALL(irep-T_thres+1,ii,:,6) = (squeeze(mean(bvardfm(:,ii,:),1))' - Y(irep+ii,:)).^2;
%        msfe_ALL(irep-T_thres+1,ii,:,7) = (squeeze(mean(bvarbnch(:,ii,:),1))' - Y(irep+ii,:)).^2;
        msfe_ALL(irep-T_thres+1,ii,:,8) = (squeeze(mean(barbnch(:,ii,:),1))' - Y(irep+ii,:)).^2;
            
       if ndraws > 1    
           for j = series_to_eval
               PL(irep-T_thres+1,ii,j,1) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),1),yi(:,j));
               PL(irep-T_thres+1,ii,j,2) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),2),yi(:,j));
               PL(irep-T_thres+1,ii,j,3) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),3),yi(:,j));
               PL(irep-T_thres+1,ii,j,4) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),4),yi(:,j));
               PL(irep-T_thres+1,ii,j,5) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),5),yi(:,j));
               PL(irep-T_thres+1,ii,j,6) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),6),yi(:,j));
               PL(irep-T_thres+1,ii,j,7) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),7),yi(:,j)); 
               PL(irep-T_thres+1,ii,j,8) = ksdensity(squeeze(fore(irep-T_thres+1,:,ii,j),8),yi(:,j)); 
           end
       else
           PL(irep-T_thres+1,ii,series_to_eval,1) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,2) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,3) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,4) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,5) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,6) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,7) = NaN;
           PL(irep-T_thres+1,ii,series_to_eval,8) = NaN;
       end
    end
end

save([pwd,'/Output/',sprintf('%s_%s_%g_%g_%g_%g_%g_%g.mat','FORECASTING',VAR_size,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp)],'Y','fore*','msfe*','PL*','-mat');