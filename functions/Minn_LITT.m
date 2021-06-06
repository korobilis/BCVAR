function [a_prior,V_prior,iV_prior,Sigma_0] = Minn_LITT(Y,Ylag,par)
%==========================================================================
% Function calculates the Litterman (Minnesota) prior
%==========================================================================
% INPUTS:
%    Y: Dependent variables
%    Ylag: lags of dependent variables (no intercept)
%    par: structure array with tunning parameters:
%                  par.gamma:  Tunning parameter for mean of own lag coeffs
%                             (usually 0 for stationary data, 1 for RW prior)
%                  par.alpha:  Tunning parameter for variance of intercepts
%                             (usually set to large number, e.g. 100)
%                  par.lambda: Crucial parameter for shrinkage in variance
%                             (usually optimize in a grid)
%                  par.theta:  Additional shrinkage parameter for coeffs of
%                              variable i appearing in eq j (Banbura et al 
%                              set this equal to 1, but usually 0<theta<1)
%==========================================================================

[t,M] = size(Y);
p = size(Ylag,2)./M;
K = M*(M*p+1);

% 1. Minnesota Mean on VAR regression coefficients
A_prior = [zeros(1,M); par.gamma*eye(M); zeros((p-1)*M,M)]';
a_prior = A_prior(:);

% 2. Minnesota Variance on VAR regression coefficients
% Now get residual variances of univariate p-lag autoregressions. Here
% we just run the AR(p) model on each equation, ignoring the constant
% and exogenous variables (if they have been specified for the original
% VAR model)
p_MIN = p;
sigma_sq = zeros(M,1); % vector to store residual variances
for i = 1:M
    Ylag_i = Ylag;
    X_i = [ones(t-p_MIN+p,1) Ylag_i(:,i:M:M*p_MIN)];
    Y_i = Y(:,i);
    % OLS estimates of i-th equation
    alpha_i = (X_i'*X_i)\(X_i'*Y_i);
    sigma_sq(i,1) = (1./(t-p_MIN+1))*(Y_i - X_i*alpha_i)'*(Y_i - X_i*alpha_i);
end

% Now define prior hyperparameters.
% Create an array of dimensions K x M, which will contain the K diagonal   
% elements of the covariance matrix, in each of the M equations.
V_i = zeros(K/M,M);
    
% index in each equation which are the own lags
ind = zeros(M,p);
for i=1:M
    ind(i,:) = i+1:M:K/M;
end
for i = 1:M  % for each i-th equation
    for j = 1:K/M   % for each j-th RHS variable   
        if j==1 % if there is constant, use this code
            V_i(j,i) = par.alpha;%*sigma_sq(i,1); % variance on intercept               
        elseif find(j==ind(i,:))>0
            V_i(j,i) = par.lambda./((ceil(j-1)/M).^2); % variance on own lags           
            % Note: the "ceil((j-1)/M)" command finds the associated lag 
            % number for each parameter
        else
            for kj=1:M
                if find(j==ind(kj,:))>0
                    ll = kj;                   
                end
            end    % variance on other lags
            V_i(j,i) = (par.lambda*par.theta*sigma_sq(i,1))./(((ceil(j-1)/M).^2)*sigma_sq(ll,1));           
        end        
    end
end

% Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
% V_i_T = V_i';
V_prior = single(diag(V_i(:)));     % this is the prior variance of the vector alpha
iV_prior = single(diag(1./V_i(:)));  

Sigma_0 = single(diag(sigma_sq));
