function [Y_forc] = BVAR_OLS(Y,p,h,ndraws)
% BAR: Bayesian Univariate AutoRegression
% Based on Koop, Korobilis and Pettenuzzo (2015)
% =========================================================================
% This code estimates the Bayesian Compressed Vector Autoregression model
% (will expand)
% INPUTS         Y: TxM matrix of data,
%                p: number of lags,
%                h: forecast horizons,
%                ndraws: number of Monte Carlo draws
%                Note: ndraws = 1 gives point forecasts
% OUTPUT: Y_forc is the predictive density
%
% =========================================================================
% Written by Davide Pettenuzzo, and Dimitris Korobilis
% Date: October 5, 2015
% =========================================================================

% ======================| Set up and estimate VAR |========================
[T, M] = size(Y);       % Dimensions of VAR

% Take lags, and correct observations
Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
T = T-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
Y = Y(p+1:end,:);

% ==============| Start computation

% Storage matrices
Y_forc = zeros(ndraws,h,M);

X = [ones(T,1),Ylag];
    
% Analytical parameter posteriors
A_OLS = (X'*X)\(X'*Y);
SIGMA_OLS = (Y - X*A_OLS)'*(Y - X*A_OLS)./T;

% Matrices in companion form
By = [A_OLS(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
Sy = zeros(M*p,M*p);
Sy(1:M,1:M) = SIGMA_OLS;
miu = zeros(M*p,1);
miu(1:M,:) = A_OLS(1,:)';
    
% Prediction
Y_pred = zeros(ndraws,h,M); % Matrix to save prediction draws
    
% Now do prediction using standard formulas (see Lutkepohl, 2005)
VAR_MEAN = 0;
VAR_VAR = 0;
X_FORE = [Y(end,:) X(end,2:M*(p-1)+1)];
BB = eye(M*p);
for ii = 1:h % not very efficient, By^(ii-1) can be defined once
    FORECASTS = VAR_MEAN + (BB*By)*X_FORE';
    if ndraws > 1
        VAR_VAR = VAR_VAR + BB*Sy*BB';
        Y_pred(:,ii,:) = (repmat(FORECASTS(1:M),1,ndraws) +  chol(VAR_VAR(1:M,1:M))'*randn(M,ndraws))';
    else
        Y_pred(:,ii,:) = repmat(FORECASTS(1:M),1,ndraws);
    end
    BB = BB*By;
end

% Store predictive draws/mean
Y_forc = Y_pred;
