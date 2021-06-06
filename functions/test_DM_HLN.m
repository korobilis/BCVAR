function [DM_adj,pval_L,pval_LR,pval_R] = test_DM_HLN(e_bench_sq,e_model_sq,h)

% @DMariano( options ) actual f1 f2 start end
% Diebold-Mariano forecast comparison test with Harvey et al (1997)
% adjustment
%
% Reference:
%   Diebold, Francis X., and Roberto S. Mariano(1995), "Comparing
%    Predictive Accuracy," JBES, v.13, no.3, pp. 253-63.
%
%   Harvey, Leybourne, Newbold(1997), "Testing the equality of prediction mean
%    squared errors", International Journal of Forecasting, vol 13, pp 289-291.
%


% Compute the gap
d = e_bench_sq-e_model_sq;
T = size(d,1);

% The Diebold-Mariano statistic is simply the "t-stat" on a regression of d on a
% constant with HAC standard errors.
var=NeweyWest(d,h);
tstats=sqrt(T)*mean(d)/sqrt(var);

% Harvey et al adjustment
HLN_adjustment = (T+1-2*h+(h*(h-1)/T))/T;
DM_adj         = tstats*sqrt(HLN_adjustment);


%one sided test H0: MSE_bench<=MSE_model, H1: MSE_bench>MSE_model
%          (model not better than bench), (model better than bench)
pval_R=1-normcdf(DM_adj,0,1);

%one sided test H0: MSE_bench>=MSE_model, H1: MSE_bench<MSE_model
%           (model not worse than bench), (model worse than bench)
pval_L=normcdf(DM_adj,0,1);

%two sided test H0: MSE_bench=MSE_model, H1:MS_bench different from MSE_model
%                   (model equal bench), (model different from bench)
% Compute pvalue only if DM_adj>0, otherwise if DM_adj<=0, set pvalue to 0
if DM_adj>0;
    pval_LR=(1-normcdf(DM_adj,0,1))+normcdf(-DM_adj,0,1);
elseif DM_adj<=0 
    %pval_LR=(1-normcdf(-DM_adj,0,1))+normcdf(DM_adj,0,1);
    pval_LR=1;
end

