function [Y_f_BMA] = BDFM(Yraw,p,constant,h,ndraws,m_max,series_to_eval)
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
% Date: 19 March, 2015
% =========================================================================

% ======================| Set up and estimate VAR |========================
[T, M] = size(Yraw);     % Dimensions of VAR
K = M*p + constant;      % Number of VAR coefficients in each equation

% Take lags, and correct observations
Ylag = mlag2(Yraw,1);
Ylag = Ylag(1+1:end,:);
T = T-1;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
X = [ones(T,1) Ylag]; X = X(:,2-constant:end);
Y = Yraw(1+1:end,:);
% % 2/ VAR matrices for SUR form
% x = kron(eye(M),X);
% y = Y(:);

% ==============| Get OLS coefficients
% beta_OLS = (x'*x)\(x'*y);
% err = reshape((y-x*beta_OLS),T,M);
% SIGMA_OLS = err'*err/(T-M);
% iSIGMA_OLS = randn(M,M);
% Winv = kron(iSIGMA_OLS,eye(T));
% index_kron = find(Winv~=0);

% ==============| Priors for FAVAR
m_max = min(m_max,size(Yraw,2)-1);
m_min = 1;                    % Minimum number of factors (just leave to 1)

% Storage matrices
BIC       = NaN(m_max-m_min+1,1);
BIC4lags  = NaN(p,1);
Y_forc    = zeros(ndraws,M,m_max-m_min+1,h);

disp(['Now running BDFM']);

nt=size(Yraw,1);
nn=size(Yraw,2);
cn=sqrt(min([nt;nn]));
clear('icp1','icp2','icp3');

for irep = m_min:m_max
%     if mod(irep,round(((m_max-m_min+1))/20))==0     
%         fprintf('%d %% \t',100*(irep/((m_max-m_min+1))));
%         %disp([num2str(100*(irep/((m_max-m_min+1)))) '% completed'])
%     end
    
    % Extract factors and create VAR matric   
    [F,L,~,~] = pc_factor(standard(Y),irep);
    YF = F;   
    YX = Y;    
    m  = irep;
    L = (olssvd(YX,[ones(T,1) YF]))';    
    e = YX - [ones(T,1) YF]*L';
    sigma2 = diag(diag(e'*e./T));
    sigma1 = sqrt(sigma2);
    beta_OLS  = cell(p,1);
    sigmaf    = cell(p,1);
    for nlag = 1:p
        K = m*nlag+constant;  

        % Obtain the errors of the VAR(1) equation
        yy = YF(nlag+1:end,:);
        xx = mlag2(YF,nlag); xx = [ones(T-nlag,1) xx(nlag+1:end,:)];
        beta_OLS{nlag,1} = (xx'*xx)\(xx'*yy);
        sigmaf{nlag,1} = (yy - xx*beta_OLS{nlag,1})'*(yy - xx*beta_OLS{nlag,1})/(T-nlag-1);
%         fac_cov = L(:,2:end)*sigmaf{nlag,1}*L(:,2:end)' + sigma2;
%         fac_pars = beta_OLS{nlag,1}*L(:,2:end)';
        BIC4lags(nlag,1) = log(det(sigmaf{nlag,1})) + (log(T)/T)*numel(beta_OLS{nlag,1});
    end
    
    % Optimize lags using BIC
    [value,best_lag] = min(BIC4lags);
    By = [beta_OLS{best_lag,1}(1+constant:end,:)'; eye(m*(best_lag-1)) zeros(m*(best_lag-1),m)];
    Sy = zeros(m*best_lag,m*best_lag);
    Sy(1:m,1:m) = sigmaf{best_lag,1};
    if constant              
        miu = zeros(m*best_lag,1);
        miu(1:m,:) = beta_OLS{best_lag,1}(1,:)';
    end
                  
    % -------| STEP 3: Prediction
    Y_pred = zeros(ndraws,M,h); % Matrix to save prediction draws
            
    % Now do prediction using standard formulas (see Lutkepohl, 2005)
    VAR_MEAN = 0;
    VAR_VAR = 0;
    X_FORE = [yy(end,:) xx(end,1+constant:m*(best_lag-1)+constant)];
    BB = speye(m*best_lag);
    for ii = 1:h % not very efficient, By^(ii-1) can be defined once
        if constant
            VAR_MEAN =  VAR_MEAN + BB*miu;
        end
        FORECASTS = VAR_MEAN + (BB*By)*X_FORE';
        if ndraws > 1
            VAR_VAR = VAR_VAR + BB*Sy*BB';
            Y_pred(:,:,ii) =  repmat(L(:,1)',ndraws,1) + ( L(:,2:end)*(repmat(FORECASTS(1:m),1,ndraws) +  chol(VAR_VAR(1:m,1:m))'*randn(m,ndraws)) )' + randn(ndraws,M)*sigma1;
        else
            Y_pred(:,:,ii) = repmat(L(:,1)',ndraws,1) + ( L(:,2:end)*repmat(FORECASTS(1:m),1,ndraws) )';
        end
        BB = BB*By;
    end
    
    % Store predictive draws/mean   
    Y_forc(:,:,irep,:) = Y_pred(:,:,:);
    % -------| STEP 4: Calculate BIC and convert these to model weights
    BIC(irep,1) = value;
end

% % -------| STEP 4: Calculate BIC and convert these to model weights
% Post_weights = exp(-.5*(BIC - min(BIC))) / sum(exp(-.5*(BIC - min(BIC)))); % - min(BIC) in both numerator and denominator is for stability
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
