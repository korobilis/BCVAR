function [Y_pred] = FAVARpredsim(Y,X,M,p,constant,L,B,SIGMA,h,ndraws)
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
J = max(size(L));
Y_pred = zeros(ndraws,J,h); % Matrix to save prediction draws

X_fore = [1 Y(end,:) X(end,1+constant:M*(p-1)+constant)]; 
X_fore = X_fore(:,2-constant:end);       
for irep = 1:ndraws               
    % Forecast of T+1 conditional on data at time T
    Y_hat = X_fore*B(:,:,irep) + randn(1,M)*chol(SIGMA);
    Y_pred(irep,:,1) = Y_hat*L';                       
       
    for i = 1:h-1  % Predict T+2, T+3 until T+h
        if i <= p                       
            X_new_temp = [1 Y_hat X_fore(:,1+constant:M*(p-i)+constant)];
            X_new_temp = X_new_temp(:,2-constant:end);
            % This gives the forecast T+i for i=1,..,p
            Y_temp = X_new_temp*B(:,:,irep) + randn(1,M)*chol(SIGMA);
            Y_pred(irep,:,i+1) = Y_temp*L';
            
            Y_hat = [Y_temp, Y_hat];
        else
            X_new_temp = [1 Y_hat(:,1:M*p)];
            X_new_temp = X_new_temp(:,2-constant:end);
            Y_temp =  X_new_temp*B(:,:,irep) + randn(1,M)*chol(SIGMA);
            Y_pred(irep,:,i+1) = Y_temp*L';
            
            Y_hat = [Y_temp, Y_hat];
        end
    end %  the last value of 'Y_temp' is the prediction T+h
end