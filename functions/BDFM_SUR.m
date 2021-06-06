function [Y_f_BMA] = BDFM_SUR(Yraw,p,constant,h,ndraws,m_max)
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
Ylag = mlag2(Yraw,p);
Ylag = Ylag(p+1:end,:);
T = T-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
X = [ones(T,1) Ylag]; X = X(:,2-constant:end);
Y = Yraw(p+1:end,:);
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

% Prior variance of beta
K_max = m_max*p + constant;
SIGMA_beta = 10*eye(K_max*m_max);      % Prior variance of VAR coefficients
i_SIGMA_beta = 1/10*eye(K_max*m_max);  % Inverse prior variance of VAR coefficients

% Storage matrices
beta_OLS  = cell(m_max-m_min+1,1);
beta_mean = cell(m_max-m_min+1,1);
beta_var  = cell(m_max-m_min+1,1);
BIC       = NaN(m_max-m_min+1,1);
Y_forc    = zeros(ndraws,M,m_max-m_min+1,h);

disp(['Now running BDFM_SUR']);

nt=size(Yraw,1);
nn=size(Yraw,2);
cn=sqrt(min([nt;nn]));
clear('icp1','icp2','icp3');

for irep = m_min:m_max
%     if mod(irep,round(((m_max-m_min+1))/20))==0     
%         fprintf('%d %% \t',100*(irep/((m_max-m_min+1))));
%         %disp([num2str(100*(irep/((m_max-m_min+1)))) '% completed'])
%     end
    
    % Extract factors and create VAR matrices
    [F,L,~,ssr] = pc_factor(Yraw,irep);
%     v=ssr/(nt*nn);
%     icp1(irep,1) = log(v)+irep*((nt+nn)/(nn*nt))*log((nn*nt)/(nn+nt));
%     icp2(irep,1) = log(v)+irep*((nt+nn)/(nn*nt))*log(cn*cn);
%     icp3(irep,1) = log(v)+irep*(1/(cn*cn))*log(cn*cn);
    
    K = irep*p+constant;
    % take lags
    Flag = mlag2(F,p);
    Flag = Flag(p+1:end,:);
    % 1/ VAR matrices for traditional matrix form
    FX = [ones(size(Flag,1),1) Flag]; FX = FX(:,2-constant:end);
    FY = F(p+1:end,:);
    % 2/ VAR matrices for SUR form
    fx = kron(eye(irep),FX);
    fy = FY(:);
    
    % OLS quantities for FAVAR
    beta_OLS{irep,1} = (fx'*fx)\(fx'*fy);
    err = reshape((fy-fx*beta_OLS{irep,1}),T,irep);
    SIGMA_OLS = err'*err/(T-irep);
    iSIGMA_OLS = inv(SIGMA_OLS);
    Winv = kron(iSIGMA_OLS,eye(T));
    index_kron = find(Winv~=0);
                       
    % -------| STEP 2: Analytical parameter posteriors conditional on prior  
    % Now get Bayesian posterior for beta
    fxW = fx'*Winv;
    beta_var{irep,1} = inv(fxW*fx + i_SIGMA_beta(1:K*irep,1:K*irep));
    beta_mean{irep,1} = beta_var{irep,1}*fxW*fy;       
    
    % -------| STEP 3: Prediction
    if ndraws ~= 1
        % ===| OPTION 1: Predictive density (simulation) |===
        beta_draw = repmat(beta_mean{irep,1},1,ndraws) + chol(beta_var{irep,1})'*randn(K*irep,ndraws); % draws from compressed VAR model
        B = reshape(beta_draw,K,irep,ndraws);  % draws from original VAR model, in matrix form
    
        [Y_pred] = FAVARpredsim(FY,FX,irep,p,constant,L,B,SIGMA_OLS,h,ndraws);
    elseif ndraws == 1
        % ===| OPTION 2: Point forecast |===
        A_mean = reshape(beta_mean{irep,1},K,irep);
        Y_pred = zeros(ndraws,M,h); % Matrix to save prediction draws

        FX_fore = [1 FY(end,:) FX(end,1+constant:irep*(p-1)+constant)]; 
        FX_fore = FX_fore(:,2-constant:end);                 
        FY_hat = FX_fore*A_mean;
        Y_pred(:,:,1) = FY_hat*L';                       
        for i = 1:h-1  % Predict T+2, T+3 until T+h       
            if i <= p                              
                FX_new_temp = [1 FY_hat FX_fore(:,1+constant:irep*(p-i)+constant)];
                FX_new_temp = FX_new_temp(:,2-constant:end);
                FY_temp = FX_new_temp*A_mean;       
                Y_pred(:,:,i+1) = FY_temp*L';            
                FY_hat = [FY_temp, FY_hat];
            else
                FX_new_temp = [1 FY_hat(:,1:irep*p)];
                FX_new_temp = FX_new_temp(:,2-constant:end);
                FY_temp =  FX_new_temp*A_mean;
                FY_pred(:,:,i+1) = FY_tempL';        
                FY_hat = [FY_temp, FY_hat];
            end
        end
    end
    
    % Store predictive draws/mean
    Y_forc(:,:,irep,:) = Y_pred;
    
    % -------| STEP 4: Calculate BIC and convert these to model weights
    err_post = reshape((fy - fx*beta_mean{irep,1}),T,irep);
    BIC(irep,1) = log(det(err_post'*err_post/T)) + (log(T)/T)*numel(beta_mean{irep,1});
end

% disp([[{'Num fact'},{'ICP1'},{'ICP2'},{'ICP3'}];num2cell([(m_min:m_max)',icp1,icp2,icp3])]);
% [ii,jj] = min([icp1 icp2 icp3]);
% disp('ICP optimal number of factors');
% disp([[{'ICP1'},{'ICP2'},{'ICP3'}];num2cell([jj])]);

% -------| STEP 4: Calculate BIC and convert these to model weights
% BIC       = Inf*ones(m_max-m_min+1,1);
% BIC(jj(3)) = 0; 
Post_weights = exp(-.5*(BIC - min(BIC))) / sum(exp(-.5*(BIC - min(BIC)))); % - min(BIC) in both numerator and denominator is for stability

% Compute BMA predictions
Y_f_BMA = zeros(ndraws,h,M);
weights = repmat(Post_weights,[1,M,ndraws]); weights = permute(weights,[3 2 1]);
for ii = 1:h
    Y_f_BMA(:,ii,:) = sum( Y_forc(:,:,:,ii).*weights,3);
end
Y_f_BMA = squeeze(Y_f_BMA);

fprintf('\n');