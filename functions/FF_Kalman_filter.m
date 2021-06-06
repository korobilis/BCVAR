function [beta_update,SIGMA_t,lambda_t,delta_t] = FF_Kalman_filter(y_t,x_t,M,K,t,lambda,delta,forgetting,beta_0_prmean,beta_0_prvar,Sigma_0)

beta_pred = zeros(K,t);
beta_update = zeros(K,t);
C_ttm1 = zeros(K,K,t);
C_tt = zeros(K,K,t);
e_t = zeros(M,t);
SIGMA_t = zeros(M,M,t);

lambda_t = zeros(t,1);
delta_t = zeros(t,1);
% set minimum values attained by lambda and delta
lmin = 0.98;
dmin = 0.94;

% 1\ Kalman filter  
for irep = 1:t    
    % Predict
    if irep==1
        beta_pred(:,irep) = beta_0_prmean;
        C_ttm1(:,:,irep) = beta_0_prvar;
    else
        beta_pred(:,irep) = beta_update(:,irep-1);
        C_ttm1(:,:,irep) = (1./lambda)*C_tt(:,:,irep-1);     
    end

    e_t(:,irep) = y_t(irep,:)' - x_t((irep-1)*M+1:irep*M,:)*beta_pred(:,irep);
   
    e2_t = e_t(:,irep)*e_t(:,irep)';
    if irep==1
        SIGMA_t(:,:,irep) = delta*Sigma_0;   
    else
        SIGMA_t(:,:,irep) = delta*SIGMA_t(:,:,irep-1) + (1-delta)*e2_t;
    end
    
    %update beta[t]
    Rx = C_ttm1(:,:,irep)*x_t((irep-1)*M+1:irep*M,:)';
    KV = SIGMA_t(:,:,irep) + x_t((irep-1)*M+1:irep*M,:)*Rx;
    KG = Rx/KV;       
    beta_update(:,irep) = beta_pred(:,irep) + KG*e_t(:,irep); %#ok<*MINV>
    C_tt(:,:,irep) = C_ttm1(:,:,irep) - KG*(x_t((irep-1)*M+1:irep*M,:)*C_ttm1(:,:,irep));
    
    % Update forgetting/decay factors
    if forgetting == 1
        % update forgetting factor
        lambda_t(irep)  = lmin + (1-lmin)*exp(-.5*abs(e2_t/KV));
        lambda = lambda_t(irep);
        
        % update decay factor
        if irep>12
            delta_t(irep) = dmin + (1-dmin)*exp(-.5*kurtosis(e_t(:,irep-12:irep)));
        else
            delta_t(irep) = delta;
        end
        delta = delta_t(irep);                
    end
end