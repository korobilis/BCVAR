function [Y_pred] = RVARpredsim(Yraw,F,M,p,constant,beta_fore,SIGMA,h,ndraws)
%==========================================================================
% Function does iterated forecasts using predictive simulation
%==========================================================================
% INPUTS:
%    Y: Dependent variables
%    X: RHS variables (mainly lags)
%    M: # of endogenous variables
%    p: # of lags
%    constant: 0 no intercept; 1 intercept
%    B: M x (M+1)p x ndraws matrix of VAR coefficients
%    SIGMA: M x M covariance matrix (OLS point estimate in our case)
%    h: Forecast horizon
%    ndraws: # of draws for predictive simulation
%==========================================================================

Y_pred = zeros(ndraws,M,h); % Matrix to save prediction draws
nfac = size(F,2);

xx = [1 F(end,:)];
xx = xx(:,2-constant:end);
x_fore = kron(eye(M),xx);
Y_hat = x_fore*beta_fore + chol(SIGMA)'*randn(M,ndraws);
Y_new = [repmat(Yraw,1,1,ndraws); reshape(Y_hat,1,M,ndraws)];
Ylg = mlag2(Y_new(:,:,1),p); Ylg = Ylg(p+1:end,:);
[fn,~,~,~] = pc_factor(standard(Ylg),nfac);  
for i = 1:h
    xx = [1 fn(end,:)];
    xx = xx(:,2-constant:end);
    x_fore = kron(eye(M),xx);
    Y_hat = x_fore*beta_fore + chol(SIGMA)'*randn(M,ndraws);
    Y_new = [repmat(Yraw,1,1,ndraws); reshape(Y_hat,1,M,ndraws)];
    Ylg = mlag2(Y_new(:,:,1),p); Ylg = Ylg(p+1:end,:);
    [fn,~,~,~] = pc_factor(standard(Ylg),nfac);
    Y_pred(:,:,i) = Y_hat';
end