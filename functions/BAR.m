function [Y_forc] = BAR(Y,p,h,ndraws)
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

% ==============| Priors for AR(p) coeffs and variance
V_beta = 10*eye(p+1);
iV_beta = 1/10*eye(p+1);  % Inverse prior variance of VAR coefficients

v_prior  = 0;
s2_prior = 1;

% ==============| Start computation

% Storage matrices
Y_forc = zeros(ndraws,h,M);


for m = 1:M
    y = Y(:,m);
    this_indx = repmat(m,1,p) + (0:M:M*(p-1)); % selects lags for m series
    X = [ones(T,1),Ylag(:,this_indx)];
    
    % Analytical parameter posteriors
    A_OLS = (X'*X)\(X'*y);
    SSE = (y - X*A_OLS)'*(y - X*A_OLS);
    V_post = inv( iV_beta + X'*X );
    A_post = V_post*(X'*X*A_OLS);
    
    S_post = SSE + v_prior*s2_prior + A_OLS'*inv(V_beta+inv(X'*X))*A_OLS;
    v_post = T + v_prior;
    
    % Matrices in companion form
    By = [A_post(2:end,:)'; eye(p-1) , zeros(p-1,1)];
    Sy = zeros(p,p);
    Sy(1,1) = S_post./v_post;
    miu = zeros(p,1);
    miu(1,:) = A_post(1,:)';
    
    
    % Prediction
    Y_pred = zeros(ndraws,h); % Matrix to save prediction draws
    
    % Now do prediction using standard formulas (see Lutkepohl, 2005)
    VAR_MEAN = 0;
    VAR_VAR = 0;
    X_FORE = [y(end,:) X(end,2:p)];
    BB = eye(p);
    for ii = 1:h % not very efficient, By^(ii-1) can be defined once
        
        VAR_MEAN =  VAR_MEAN + BB*miu;
        
        FORECASTS = VAR_MEAN + (BB*By)*X_FORE';
        if ndraws > 1
            VAR_VAR = VAR_VAR + BB*Sy*BB';
            Y_pred(:,ii) = (repmat(FORECASTS(1),1,ndraws) +  chol(VAR_VAR(1,1))'*randn(1,ndraws))';
        else
            Y_pred(:,ii) = repmat(FORECASTS(1),1,ndraws);
        end
            
        BB = BB*By;
    end
    
    
    % Store predictive draws/mean
    Y_forc(:,:,m) = Y_pred;

end