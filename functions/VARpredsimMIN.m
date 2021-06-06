function [Y_pred] = VARpredsimMIN(Y,X,M,p,constant,B,SIGMA,h,ndraws)
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

% -------| Initialize matrix to store predicted density
Y_pred = zeros(ndraws,M,h); % Matrix to save prediction draws
S = chol(SIGMA);

X_fore = [Y(end,:) X(end,1+constant:M*(p-1)+constant) 1]; 
X_fore = X_fore(:,1:end--(constant==0));       
for irep = 1:ndraws               
    % Forecast of T+1 conditional on data at time T
    Y_hat = X_fore*B(:,:,irep) + randn(1,M)*S;
    Y_pred(irep,:,1) = Y_hat;                       
       
    for i = 1:h-1  % Predict T+2, T+3 until T+h
        if i <= p                       
            X_new_temp = [Y_hat X_fore(:,1+constant:M*(p-i)+constant) 1];
            X_new_temp = X_new_temp(:,1:end-(constant==0));
            % This gives the forecast T+i for i=1,..,p
            Y_temp = X_new_temp*B(:,:,irep) + randn(1,M)*S;
            Y_pred(irep,:,i+1) = Y_temp;
            
            Y_hat = [Y_temp, Y_hat];
        else
            X_new_temp = [Y_hat(:,1:M*p) 1];
            X_new_temp = X_new_temp(:,1:end-(constant==0));
            Y_temp =  X_new_temp*B(:,:,irep) + randn(1,M)*S;
            Y_pred(irep,:,i+1) = Y_temp;
            
            Y_hat = [Y_temp, Y_hat];
        end
    end %  the last value of 'Y_temp' is the prediction T+h
end