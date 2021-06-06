function [MHLN,pval_L,pval_LR,pval_R]=hlntest_DP(e_bench, e_model,h)

%e_bench=y-y_bench;
%e_model=y-y_model;

d=(e_bench-e_model).*e_bench;
dhat=mean(d);

% Calculate the variance of the loss differential, taking into account autocorrelation.
T = size(e_bench,1);
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


phi=sum((d-dhat).^2)/length(d);
Vdhat=phi/length(d);

keyboard
[tstat,beta,Vdhat_NW] = nwest(d,ones(length(d),1),h-1);

MHLN=((length(d)-1)/(length(d)))*dhat/sqrt(Vdhat);
%pval=1-tcdf(MHLN,length(d)-1); %one side test

%one sided test H0: MSE1=MSE2, H1=MSE1>MSE2
pval_R=1-tcdf(MHLN,length(d)-1);
%one sided test H0: MSE1=MSE2, H1=MSE1<MSE2
pval_L=tcdf(MHLN,length(d)-1);
%two side test H0: MSE1=MSE2, H1=MS1 different from MSE2
if MHLN>0;
    pval_LR=(1-tcdf(MHLN,length(d)-1))+tcdf(-MHLN,length(d)-1);
elseif MHLN<=0 || isnan(MHLN)
    pval_LR=(1-tcdf(-MHLN,length(d)-1))+tcdf(MHLN,length(d)-1);
end





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

function [tstat,beta,nwerr]=nwest(y,x,nlag)

% PURPOSE: computes Newey-West adjusted heteroscedastic-serial
%          consistent Least-squares Regression
%---------------------------------------------------
% USAGE: results = nwest(y,x,nlag)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar)
%     nlag = lag length to use
% References:  Gallant, R. (1987),
%  "Nonlinear Statistical Models," pp.137-139.
%---------------------------------------------------

% written by James P. LeSage and modified by Lutz Kilian
[nobs nvar] = size(x);
xpxi = inv(x'*x);
beta    = xpxi*(x'*y);
yhat    = x*beta;
resid   = y - yhat;
sigu = resid'*resid;
sige    = sigu/(nobs-nvar);

% perform Newey-West correction
emat = [];
for i=1:nvar;
    emat = [emat
        resid'];
end;

hhat=emat.*x';
G=zeros(nvar,nvar); w=zeros(2*nlag+1,1);
a=0;

while a~=nlag+1;
    ga=zeros(nvar,nvar);
    w(nlag+1+a,1)=(nlag+1-a)/(nlag+1);
    za=hhat(:,(a+1):nobs)*hhat(:,1:nobs-a)';
    if a==0;
        ga=ga+za;
    else
        ga=ga+za+za';
    end;
    G=G+w(nlag+1+a,1)*ga;
    a=a+1;
end;

V=xpxi*G*xpxi;
nwerr= sqrt(diag(V));

tstat = beta./nwerr; % Newey-West t-statistics
