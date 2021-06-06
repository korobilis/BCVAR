function [Y_f_BMA] = BCTRVAR_CONJ(Y,p,constant,h,ndraws,RP_type,n_psi,stdata,apply_bcr,weight_scheme,cov_comp,sparsity,series_to_eval)                   
% BCTRVAR: Bayesian Compressed TRiangular Vector AutoRegression
% Based on Koop, Korobilis and Pettenuzzo (2015)
% =========================================================================
% This code estimates the Bayesian Compressed Vector Autoregression model
% (will expand)
% INPUTS         Y: TxM matrix of data, 
%                p: number of lags,
%         constant: 0: no intercept in the VAR; 1: intercept in the VAR
%                h: forecast horizons, 
%           ndraws: number of Monte Carlo draws
%                   Note: ndraws = 1 gives point forecasts
%          RP_type: type of Random Projection matrix to generate (choice 1-3)
%            n_psi: choose how many psi's to generate (Note: only works for RP_type=1);
%           stdata: 0: do nothing; 1: standardize data;
%        apply_bcr: 1:everywhere; 2:intercepts excluded; 3: intercepts and first own lags
%    weight_scheme: 1: 1/N; 2: BIC for whole VAR; 3: equation-by-equation BIC
%         cov_comp: 0: Don't compress covariance matrix; 1: Compress cov matrix
%  
% OUTPUT: Y_forc_mean is the posterior mean of the BMA predictive density
%
% =========================================================================
% Written by Davide Pettenuzzo, and Dimitris Korobilis
% Date: 9 September, 2015
% =========================================================================
rng('shuffle');
% ======================| Set up and estimate VAR |========================
[T, M] = size(Y);       % Dimensions of VAR
N =  M*p+constant;      % Number of VAR coefficients in each M equations
K = M*N;                % Number of VAR coefficients
numa = M*(M-1)/2;       % Number of VAR covariance matrix elements

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
X = [ones(T,1) Ylag];
X = X(:,2-constant:end);
Y = Y(p+1:end,:);

% Create matrix of elements to compress
X_comp = cell(M,1);
Xu = cell(M,1);
for irep = 1:M
    Xu{irep} = [];
    if apply_bcr == 1
        X2c = X;
    elseif apply_bcr == 2   
        X2c = X(:,1+constant:end);
        if constant ~= 0
           Xu{irep} = X(:,1); 
        end
    elseif apply_bcr == 3
        X2c = X(:,1+constant:end);
        Xu{irep} = [X(:,1), X2c(:,irep)];
        X2c(:,irep) = [];
    end
    % Create now the matrix we will compres, including covariance matrix
    % (as contemporaneous regression effects)
    if cov_comp == 1
        if irep==1
            X_comp{irep} = X2c;
        else
            X_comp{irep} = [X2c, -Y(:,1:irep-1)];
        end
    elseif cov_comp == 0
        X_comp{irep} = X2c;
    end
end

% This is just to get index_cov which helps me construct the covariance
% matrix later using the vector of coefficient draws
COVM = randn(M,M);
COVM = COVM'*COVM;
COVM = chol(COVM)';
COVM = COVM/diag(diag(COVM));
index_cov = find(~tril(COVM));

% ==============| Priors for BCR
% 1. Model size (max is n_m1, min is n_m2)
n_m1 = 5*round(log(N));%min(N,T);
n_m2 = 1; %5*round(log(N));

% 2. Prior variance on VAR coefficients
prior_c = 0.5;
V_beta = prior_c*eye(n_m1+M+1);
iV_beta = (1/prior_c)*eye(n_m1+M+1);  % Inverse prior variance of VAR coefficients

% ==============| Start computation
m_max = n_m1;                       % This is the max m. Note that m_max<<p
m_min = n_m2;                       % This is the min m
Model_no = 0;                       % Just to count the model number

% Storage matrices
A_post = cell(M,1);
B_post = zeros(N,M);
C_post = zeros(numa,1);
PHI = cell(M,1);
sigma = zeros(M,1);
Y_forc = zeros(ndraws,M,n_psi,h);
BIC = zeros(n_psi,1);
BIC2 = zeros(n_psi,3);

for psi_l = 1:n_psi
    Model_no = Model_no + 1;
    if mod(Model_no,round(n_psi)/10)==0        
        fprintf('%d %% \t',100*(Model_no/n_psi));
        %disp([num2str(100*(Model_no/n_psi)) '% completed'])       
    end

    ind = 1;
    for irep = 1:M
        % -------| STEP 2: Generate PHI
        if sparsity == 1 && apply_bcr == 2 && cov_comp == 0
            m_max_s = 2;
            m_min_s = 2;
            size_PHIi = randi([m_min_s,m_max_s]);
            PHI{irep} = gensparseRP(p,size_PHIi,size(X_comp{irep},2));
        elseif sparsity == 0
            size_PHIi = randi([m_min,m_max]);
            PHI{irep} = genRP(RP_type,size_PHIi,size(X_comp{irep},2));
        else
            error('Sparse RP only available with apply_bcr=2')
        end
                
        % -------| STEP 2: Analytical parameter posteriors conditional on PHI
        if cov_comp == 1       
            X_new = [Xu{irep}'; PHI{irep}*X_comp{irep}'];    % Matrix X in eq irep. Xu contains unrestricted (not compressed) elements
        elseif cov_comp == 0
            X_new = [Xu{irep}'; PHI{irep}*X_comp{irep}'; -Y(:,1:irep-1)'];
        end
        n_vars(irep,1) = size(X_new,1);                          % Number of variables in eq irep. This changes in each eq due to the triangular structure
        PX = X_new; XP = PX'; PXXP = PX*XP;              % Define these here once, we need them several times        
        
        A_OLS = (PXXP)\(PX*Y(:,irep));                   % OLS of the compressed VAR  
        iV_post = ( iV_beta(1:n_vars(irep,1),1:n_vars(irep,1)) + PXXP ); % Posterior variance of coeffs in the compressed VAR
        A_post{irep,1} = iV_post\(PXXP*A_OLS);           % Posterior mean of coeffs in the compressed VAR

        SSE = (Y(:,irep) - XP*A_OLS)'*(Y(:,irep) - XP*A_OLS); % Sum of Squared Errors based on posterior mean of coeffs
        v_post = T;
        s_post = (SSE + (A_post{irep,1} - A_OLS)'*(PXXP/iV_post)*iV_beta(1:n_vars(irep,1),1:n_vars(irep,1))*(A_post{irep,1} - A_OLS))./v_post;
        sigma(irep,1) = s_post;                                      % Mean estimate of sigma in eq irep
        
        % equation-specific BIC
        if irep < max(series_to_eval)+1
            BIC2(Model_no,irep) = size(X_comp{irep},2)*log(T) + T*log(SSE/T);        
        end
        
        % Construct full VAR matrix
        if cov_comp == 1 
            if apply_bcr == 1  % everywhere
                Btemp = PHI{irep}'*A_post{irep};
            elseif apply_bcr == 2  % exclude intercepts
                Btemp = [A_post{irep}(1); PHI{irep}'*A_post{irep}(2:end)];
            elseif apply_bcr == 3  % exclude intercepts and first own lags           
                Btemp = [A_post{irep}(1:2); PHI{irep}'*A_post{irep}(3:end)];
                BB = Btemp(2);
                Btemp(2) = [];
                Btemp = [Btemp(1:irep); BB; Btemp(irep+1:end)];               
            end
        elseif cov_comp == 0
            if apply_bcr == 1  % everywhere
                Btemp = [PHI{irep}'*A_post{irep}(1:end-irep+1); A_post{irep}(end-irep+2:end)];           
            elseif apply_bcr == 2  % exclude intercepts
                Btemp = [A_post{irep}(1); PHI{irep}'*A_post{irep}(2:end-irep+1) ; A_post{irep}(end-irep+2:end)];
            elseif apply_bcr == 3  % exclude intercepts and first own lags           
                Btemp = [A_post{irep}(1:2); PHI{irep}'*A_post{irep}(3:end-irep+1) ; A_post{irep}(end-irep+2:end)];
                BB = Btemp(2);
                Btemp(2) = [];
                Btemp = [Btemp(1:irep); BB; Btemp(irep+1:end)];               
            end
        end
        B_post(:,irep) = Btemp(1:N);

        if irep>1 % Construct matrix of covariances from vector Btemp
           C_post(((irep-1)+(irep-3)*(irep-2)/2):ind,1) = Btemp(N+1:end);
           ind = ind + irep;
        end
    end
    
    % Reduced-form VAR matrices
    beta_full = B_post;
    Acov = eye(M); Acov(index_cov) = C_post; Acov = Acov';       
    beta_red_VAR = beta_full/Acov';              % Coeffs are here
    sigma_red_VAR = (Acov\diag(sigma))/Acov';    % Full cov matrix is here
    
    % Matrices in companion form   
    By = [beta_red_VAR(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
    Sy = zeros(M*p,M*p);
    Sy(1:M,1:M) = sigma_red_VAR;
    miu = zeros(M*p,1);
    miu(1:M,:) = beta_red_VAR(1,:)';

    while max(abs(eig(By))) > 1
        ind = 1;
        for irep = 1:M
            % -------| STEP 2: Generate PHI
            if sparsity == 1 && apply_bcr == 2 && cov_comp == 0
                m_max_s = 2;
                m_min_s = 2;
                size_PHIi = randi([m_min_s,m_max_s]);
                PHI{irep} = gensparseRP(p,size_PHIi,size(X_comp{irep},2));
            elseif sparsity == 0           
                size_PHIi = randi([m_min,m_max]);
                PHI{irep} = genRP(RP_type,size_PHIi,size(X_comp{irep},2));
            else                
                error('Sparse RP only available with apply_bcr=2')
            end
                
            % -------| STEP 2: Analytical parameter posteriors conditional on PHI    
            if cov_comp == 1       
                X_new = [Xu{irep}'; PHI{irep}*X_comp{irep}'];    % Matrix X in eq irep. Xu contains unrestricted (not compressed) elements
            elseif cov_comp == 0           
                X_new = [Xu{irep}'; PHI{irep}*X_comp{irep}'; -Y(:,1:irep-1)'];
            end
            n_vars(irep,1) = size(X_new,1);                          % Number of variables in eq irep. This changes in each eq due to the triangular structure
            PX = X_new; XP = PX'; PXXP = PX*XP;              % Define these here once, we need them several times        
            A_OLS = (PXXP)\(PX*Y(:,irep));                   % OLS of the compressed VAR  
            iV_post = ( iV_beta(1:n_vars(irep,1),1:n_vars(irep,1)) + PXXP ); % Posterior variance of coeffs in the compressed VAR
            A_post{irep,1} = iV_post\(PXXP*A_OLS);           % Posterior mean of coeffs in the compressed VAR
            
            SSE = (Y(:,irep) - XP*A_OLS)'*(Y(:,irep) - XP*A_OLS); % Sum of Squared Errors based on posterior mean of coeffs
            v_post = T;
            s_post = (SSE + (A_post{irep,1} - A_OLS)'*(PXXP/iV_post)*iV_beta(1:n_vars(irep,1),1:n_vars(irep,1))*(A_post{irep,1} - A_OLS))./v_post;
            sigma(irep,1) = s_post;                                      % Mean estimate of sigma in eq irep
        
            % equation-specific BIC
            BIC2(Model_no,irep) = size(X_comp{irep},2)*log(T) + T*log(SSE/T);        
                   
            % Construct full VAR matrix
            if cov_comp == 1            
                if apply_bcr == 1  % everywhere
                    Btemp = PHI{irep}'*A_post{irep};
                elseif apply_bcr == 2  % exclude intercepts               
                    Btemp = [A_post{irep}(1); PHI{irep}'*A_post{irep}(2:end)];
                elseif apply_bcr == 3  % exclude intercepts and first own lags           
                    Btemp = [A_post{irep}(1:2); PHI{irep}'*A_post{irep}(3:end)];
                    BB = Btemp(2);
                    Btemp(2) = [];
                    Btemp = [Btemp(1:irep); BB; Btemp(irep+1:end)];               
                end
            elseif cov_comp == 0           
                if apply_bcr == 1  % everywhere
                    Btemp = [PHI{irep}'*A_post{irep}(1:end-irep+1); A_post{irep}(end-irep+2:end)];           
                elseif apply_bcr == 2  % exclude intercepts              
                    Btemp = [A_post{irep}(1); PHI{irep}'*A_post{irep}(2:end-irep+1) ; A_post{irep}(end-irep+2:end)];
                elseif apply_bcr == 3  % exclude intercepts and first own lags                
                    Btemp = [A_post{irep}(1:2); PHI{irep}'*A_post{irep}(3:end-irep+1) ; A_post{irep}(end-irep+2:end)];
                    BB = Btemp(2);
                    Btemp(2) = [];
                    Btemp = [Btemp(1:irep); BB; Btemp(irep+1:end)];               
                end
            end
            B_post(:,irep) = Btemp(1:N);

            if irep>1 % Construct matrix of covariances from vector Btemp
                C_post(((irep-1)+(irep-3)*(irep-2)/2):ind,1) = Btemp(N+1:end);
                ind = ind + irep;
            end
        end
        
        % Reduced-form VAR matrices
        beta_full = B_post;
        Acov = eye(M); Acov(index_cov) = C_post; Acov = Acov';       
        beta_red_VAR = beta_full/Acov';              % Coeffs are here
        sigma_red_VAR = (Acov\diag(sigma))/Acov';    % Full cov matrix is here           
        
        % Matrices in companion form   
        By = [beta_red_VAR(2:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];
        Sy = zeros(M*p,M*p);
        Sy(1:M,1:M) = sigma_red_VAR;
        miu = zeros(M*p,1);
        miu(1:M,:) = beta_red_VAR(1,:)';        
    end
    
    % -------| STEP 3: Prediction
    Y_pred = zeros(ndraws,M,h); % Matrix to save prediction draws
    
    % Now do prediction using standard formulas (see Lutkepohl, 2005)
    VAR_MEAN = 0;
    VAR_VAR = 0;
    X_FORE = [Y(end,:) X(end,1+constant:M*(p-1)+constant)];
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
    Y_forc(:,:,Model_no,:) = repmat(mm,ndraws,1,h) + repmat(ss,ndraws,1,h).*Y_pred;
      
    % -------| STEP 4: Calculate model diagnostics (BIC, ML)
    err_post = (Y - X*beta_red_VAR);
    BIC(Model_no,1) = log(det(err_post'*err_post/T)) + (log(T)/T)*sum(n_vars);  
end

% Calculate weights
if weight_scheme == 1 % 1/N      
    Post_weights = (1./n_psi)*ones(n_psi,1);
elseif weight_scheme == 2 % BIC
    PSI = BIC - min(BIC);
    Post_weights = exp(-.5*PSI) / sum(exp(-.5*PSI));
elseif weight_scheme == 3 % BIC2
    BIC = mean(BIC2,2);
    PSI = BIC - min(BIC);
    Post_weights = exp(-.5*PSI) / sum(exp(-.5*PSI));
end

% % Compute BMS predictions
% [~,jj] = max(Post_weights);
% tmp = squeeze(Y_forc(:,:,jj,:));
% Y_f_BMA = permute(tmp,[1 3 2]);

% Compute BMA predictions
Y_f_BMA = zeros(ndraws,h,M);
if ndraws > 1
    R = mnrnd(1,Post_weights,ndraws); % random draws from multinomial distribution
    [rr,cc] = find(R==1);
    for m=1:(m_max-m_min+1)*n_psi
        this_R = rr(cc==m);
        if ~isempty(this_R)
            for ii=1:h
                Y_f_BMA(this_R,ii,:) = Y_forc(this_R,:,m,ii);
            end
        end
    end
   
else
    weights = repmat(Post_weights,[1,M,ndraws]); weights = permute(weights,[3 2 1]);
    for ii = 1:h
        Y_f_BMA(:,ii,:) = sum( Y_forc(:,:,:,ii).*weights,3);
    end
end

fprintf('\n');