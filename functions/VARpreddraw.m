function [Y_pred] = VARpreddraw(Y,X,M,p,constant,B,SIGMA,h)
%==========================================================================
% Function does iterated forecasts, one draw only (for code which estimates
% coefficients with Monte Carlo integration
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

Y_pred = zeros(M,h); % Matrix to save prediction draws
S = chol(SIGMA);

X_fore = [1 Y(end,:) X(end,1+constant:M*(p-1)+constant)]; 
X_fore = X_fore(:,2-constant:end);                 
 
% Forecast of T+1 conditional on data at time T
Y_hat = X_fore*B + randn(1,M)*S;
Y_pred(:,1) = Y_hat;                       
for i = 1:h-1  % Predict T+2, T+3 until T+h       
    if i <= p                       
        X_new_temp = [1 Y_hat X_fore(:,1+constant:M*(p-i)+constant)];
        X_new_temp = X_new_temp(:,2-constant:end);
        % This gives the forecast T+i for i=1,..,p
        Y_temp = X_new_temp*B + randn(1,M)*S;
        Y_pred(:,i+1) = Y_temp;
            
        Y_hat = [Y_temp, Y_hat];
    else
        X_new_temp = [1 Y_hat(:,1:M*p)];
        X_new_temp = X_new_temp(:,2-constant:end);
        Y_temp =  X_new_temp*B + randn(1,M)*S;
        Y_pred(:,i+1) = Y_temp;
        
        Y_hat = [Y_temp, Y_hat];
    end
end