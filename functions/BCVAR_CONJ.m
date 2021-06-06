function [Y_f_BMA] = BCVAR_CONJ(Y,p,constant,h,ndraws,RP_type,n_psi,stdata,series_to_eval)
% BCVAR: Bayesian Compressed Vector AutoRegression
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
%           stdata: 0: do nothing; 1: standardize RHS data;
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

% ==============| Priors for BCR
n_m1 = 5*round(log(N));
n_m2 = 1; %1*round(log(N));

V_beta = 10*eye(n_m1+2);
iV_beta = 1/10*eye(n_m1+2);  % Inverse prior variance of VAR coefficients
   
v_prior = M;
S_prior = eye(M);
inv_S_prior = inv(S_prior);

% ==============| Start computation
m_max = n_m1;                       % This is the max m. Note that m_max<<p
m_min = n_m2;                       % This is the min m
Model_no = 0;                       % Just to count the model number

% Storage matrices
Y_forc = zeros(ndraws,M,(m_max-m_min+1)*n_psi,h);
BIC = zeros((m_max-m_min+1)*n_psi,1);

for m_l = m_min:m_max
    for psi_l = 1:n_psi        
        Model_no = Model_no + 1;        
        if mod(Model_no,round(((m_max-m_min+1)*n_psi)/10))==0  
            fprintf('%d %% \t',round(100*(Model_no/round((m_max-m_min+1)*n_psi))));
            %disp([num2str(100*(Model_no/((m_max-m_min+1)*n_psi))) '% completed'])
        end
                
        [PHI] = genRP(RP_type,m_l,N);
        
        % -------| STEP 2: Analytical parameter posteriors conditional on PHI
        PX = PHI*X'; XP = PX';       % Define these here once, we need them several times
        PXXP = PX*XP;
        A_OLS = (PXXP)\(PX*Y);
        SSE = (Y - XP*A_OLS)'*(Y - XP*A_OLS);
        V_post = inv( iV_beta(1:m_l,1:m_l) + PXXP );
        A_post = V_post*(PXXP*A_OLS);
           
        S_post = SSE + S_prior + A_OLS'*PXXP*A_OLS - A_post'*(iV_beta(1:m_l,1:m_l) + PXXP)*A_post;
        S_post = symmetric(S_post);
        v_post = T + v_prior;

        % Get coefficients in un-compressed VAR
        B_post = PHI'*A_post;
        
        % Matrices in companion form
        By = [B_post(1+constant:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];   
        Sy = zeros(M*p,M*p);
        Sy(1:M,1:M) = S_post./v_post;
        if constant       
            miu = zeros(M*p,1);
            miu(1:M,:) = B_post(1,:)';
        end
        
        while max(abs(eig(By))) > 1
            [PHI] = genRP(RP_type,m_l,N);
        
            % -------| STEP 2: Analytical parameter posteriors conditional on PHI
            PX = PHI*X'; XP = PX';       % Define these here once, we need them several times       
            PXXP = PX*XP;
            A_OLS = (PXXP)\(PX*Y);
            SSE = (Y - XP*A_OLS)'*(Y - XP*A_OLS);
            V_post = inv( iV_beta(1:m_l,1:m_l) + PXXP );
            A_post = V_post*(PXXP*A_OLS);
           
            S_post = SSE + S_prior + A_OLS'*PXXP*A_OLS - A_post'*(iV_beta(1:m_l,1:m_l) + PXXP)*A_post;
            S_post = symmetric(S_post);
            v_post = T + v_prior;

            % Get coefficients in un-compressed VAR
            B_post = PHI'*A_post;
        
            % Matrices in companion form
            By = [B_post(1+constant:end,:)'; eye(M*(p-1)) , zeros(M*(p-1),M)];   
            Sy = zeros(M*p,M*p);
            Sy(1:M,1:M) = S_post./v_post;
            if constant       
                miu = zeros(M*p,1);
                miu(1:M,:) = B_post(1,:)';
            end 
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
        
        % -------| STEP 4: Calculate BIC and convert these to model weights
        err_post = (Y - XP*A_post);
        BIC(Model_no,1) = log(det(err_post(:,series_to_eval)'*err_post(:,series_to_eval)/T)) + (log(T)/T)*numel(A_post(:,series_to_eval));
    end
    
end
PSI = BIC - min(BIC);
Post_weights = exp(-.5*PSI) / sum(exp(-.5*PSI)); % - min(BIC) in both numerator and denominator is for stability

% Compute BMA predictions
Y_f_BMA = zeros(ndraws,h,M);
if ndraws > 1
    R = mnrnd(1,Post_weights,ndraws); % random draws from multinomial distribution
    [rr,cc] = find(R==1);
    for m=1:(m_max-m_min+1)*n_psi
        this_R = rr(cc==m);
        for ii=1:h
            Y_f_BMA(this_R,ii,:) = Y_forc(this_R,:,m,ii);
        end
    end
   
else
    weights = repmat(Post_weights,[1,M,ndraws]); weights = permute(weights,[3 2 1]);
    for ii = 1:h
        Y_f_BMA(:,ii,:) = sum( Y_forc(:,:,:,ii).*weights,3);
    end
end

fprintf('\n');