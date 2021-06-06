
function [DM,pval_L,pval_LR,pval_R]= dmtest(e_bench_sq,e_model_sq,h)



% Numerator
d = (e_bench_sq) - (e_model_sq);
dMean = mean(d);

% (Denum) Ralculate the variance of the loss differential, taking into account autocorrelation.
T = size(e_bench_sq,1);
gamma0 = var(d);
if h > 1
    gamma = zeros(h-1,1);
    for i = 1:h-1
        sampleCov = cov(d(1+i:T),d(1:T-i));
        gamma(i) = sampleCov(2);
    end
    varD = gamma0 + 2*sum(gamma);
else
    varD = gamma0;
end

% Diebold-Mariano statistic DM ~N(0,1)
DM = dMean / sqrt ( (1/T)*varD );  %equivalent to R=OLS(ones(T,1),d); R.tstat==DM  
    
%one sided test H0: MSE_bench<=MSE_model, H1: MSE_bench>MSE_model
%          (model not better than bench), (model better than bench)
pval_R=1-normcdf(DM,0,1); 

%one sided test H0: MSE_bench>=MSE_model, H1: MSE_bench<MSE_model
%           (model not worse than bench), (model worse than bench)
pval_L=normcdf(DM,0,1); 

%two sided test H0: MSE_bench=MSE_model, H1:MS_bench different from MSE_model
%                   (model equal bench), (model different from bench)
if DM>0;
pval_LR=(1-normcdf(DM,0,1))+normcdf(-DM,0,1); 
elseif DM<=0 || isnan(DM)
pval_LR=(1-normcdf(-DM,0,1))+normcdf(DM,0,1);     
end
%pval_LR=[1-normcdf(abs(DM),0,1)]*2;