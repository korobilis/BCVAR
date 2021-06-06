function [Y_f_BMA] = BVAR_MINN(Y,p,constant,h,ndraws,lgrid,stdata)
% Based on Koop, Korobilis and Pettenuzzo (2015)
% =========================================================================
% % p is the number of lags, H is the max forecast horizon (forecasts are
% obtained using the direct approach), Nfac_max is the max number of
% factors considered, prior is an indicator for the priors we use (1 is for
% uninformative, 2 is for Minnesota, 3 is for natural conjugate). In all
% three cases analytical solutions are available so no need to use MCMC to
% compute posterior quantities
% NOTE: There is no constant in these regressions (a constant is included
% in the notes I wrote), but could be easily added
%
% =========================================================================
% Date: 9 September, 2015
% =========================================================================

% ======================| Set up and estimate VAR |========================
[T, M] = size(Y);     % Dimensions of VAR
K = M*p+constant;     % Number of VAR coefficients in each equation

if stdata == 1
    [Y,mm,ss] = standard(Y);
else
    mm = zeros(1,M);
    ss = ones(1,M);
end

% Take lags, and correct observations
Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
T = T-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
X = [ones(T,1) Ylag]; X = X(:,2-constant:end);
Y = Y(p+1:end,:);
% % 2/ VAR matrices for SUR form
% x = kron(eye(M),X);
% y = Y(:);

% ==============| Get OLS coefficients
% warning off
% beta_OLS = (X'*X)\(X'*Y);
% err = (Y-X*beta_OLS);
% SIGMA_OLS = err'*err/(T-M);
% iSIGMA_OLS = inv(SIGMA_OLS);
% Winv = kron(iSIGMA_OLS,eye(T));
% index_kron = find(Winv~=0);
% warning on
% ==============| Priors for MINN
% Grid for Minnesota Prior
% lgrid = [1e-5, 1e-4, 1e-3, 1e-2:1e-2:1];
n_grid = length(lgrid);
mu = 0;             % 0: No sum of coefficients prior
iRW = zeros(1,M);   % prior means (0 for stationary, 1 for RW variables)    


% Storage matrices
beta_mean = cell(n_grid,1);
BIC       = NaN(n_grid,1);
Y_forc    = zeros(ndraws,M,n_grid,h);

for m_l = 1:n_grid    
    if mod(m_l,round(n_grid/20))==0   
        fprintf('%d %% \t',round(100*(m_l/n_grid)));
        %disp([num2str(100*(m_l/n_grid)) '% completed'])
    end
        
    % Obtain Minnesota prior for model m_l
    tau = lgrid(m_l);   % shrinkage coefficient
%     [Y_pred,beta] = bvarLitt2(Y,p,tau,mu,h,iRW);
%     beta_mean{m_l,1} = beta;
    [Yd,Xd,Td] = bvarLitt(Y,p,tau,mu,iRW);
    
    % -------| STEP 2: Analytical results for posterior mean of beta and SIGMA
    ZZinv = inv(Xd'*Xd + X'*X);
    ZY    = Xd'*Yd + X'*Y;
    
    beta_mean{m_l,1} = ZZinv*ZY;
    Stilde = (Y - X*beta_mean{m_l,1})'*(Y - X*beta_mean{m_l,1});

    % -------| STEP 3: Prediction
    Y_pred = zeros(ndraws,M,h); % Matrix to save prediction draws
    
    % Matrices in companion form
    By = [beta_mean{m_l,1}(1+constant:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];   
    Sy = zeros(M*p,M*p);
    Sy(1:M,1:M) = Stilde./T;
    if constant       
        miu = zeros(M*p,1);
        miu(1:M,:) = beta_mean{m_l,1}(1,:)';
    end
    
    % Now do prediction using standard formulas (see Lutkepohl, 2005)
    VAR_MEAN = 0;
    VAR_VAR = 0;
    X_FORE = [Y(end,:) X(end,1+constant:M*(p-1)+constant)];
    X_FORE = X_FORE(:,2-constant:end);
    BB = speye(M*p);
    for ii = 1:h % not very efficient, By^(ii-1) can be defined once
        if constant
            VAR_MEAN =  VAR_MEAN + BB*miu;
        end
        FORECASTS = VAR_MEAN + (BB*By)*X_FORE';
        if ndraws > 1
            VAR_VAR = VAR_VAR + BB*Sy*BB';
            Y_pred(:,:,ii) = (repmat(FORECASTS(1:M),1,ndraws) +  chol(VAR_VAR(1:M,1:M))'*randn(M,ndraws))';
        else
            Y_pred(:,:,ii) = repmat(FORECASTS(1:M),1,ndraws);
        end
        BB = BB*By;
    end
        
    % Store predictive draws/mean
    Y_forc(:,:,m_l,:) = repmat(mm,ndraws,1,h) + repmat(ss,ndraws,1,h).*Y_pred;
    
    % -------| STEP 4: Calculate BIC and convert these to model weights
    BIC(m_l,1) = log(det(Stilde/T)) + (log(T)/T)*numel(beta_mean{m_l,1});
end

% PSI = BIC - min(BIC);
% Post_weights = exp(-.5*PSI) / sum(exp(-.5*PSI)); % - min(BIC) in both numerator and denominator is for stability
% 
% % Compute BMA predictions
% Y_f_BMA = zeros(ndraws,h,M);
% weights = repmat(Post_weights,[1,M,ndraws]); weights = permute(weights,[3 2 1]);
% for ii = 1:h
%     Y_f_BMA(:,ii,:) = sum( Y_forc(:,:,:,ii).*weights,3);
% end
% Y_f_BMA = squeeze(Y_f_BMA);
[g1,g2]=(min(BIC));
Y_best = squeeze(Y_forc(:,:,g2,:));
if ndraws > 1
    Y_f_BMA = permute(Y_best,[1 3 2]);
else
    Y_f_BMA = permute(Y_best,[3 2 1]);
end
fprintf('\n');