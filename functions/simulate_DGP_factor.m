function [dates_new,Y_sim,lam_0,lam_1] = simulate_DGP_factor(K,dates,X,SIGMA,T,p,t_beg,t_end)

% Take lags, and correct observations (dropping first p obs)
Y    = X;
Ylag = mlag2(Y,p);
Ylag(1:p,:) = NaN;

% Adjust to match the beginning and end of estimation sample
Y_new     = Y(t_beg:t_end,:);
Ylag_new  = Ylag(t_beg:t_end,:);
dates_new = dates(t_beg:t_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove obs with nans and perform estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nan_indx = sum(isnan([Y_new Ylag_new]),2)>0;
disp('Small VAR has missing values for the following dates - drop them before estimation:');
disp(datestr(dates_new(nan_indx==1)));

Y_nanfree     = Y_new(nan_indx==0,:);
Ylag_nanfree  = Ylag_new(nan_indx==0,:);
dates_nanfree = dates_new(nan_indx==0);

X_nanfree = [ones(size(Ylag_nanfree,1),1),Ylag_nanfree];
A_OLS = (X_nanfree'*X_nanfree)\(X_nanfree'*Y_nanfree); % This is the matrix of regression coefficients
SSE = (Y_nanfree - X_nanfree*A_OLS)'*(Y_nanfree - X_nanfree*A_OLS);
SIGMA_OLS = SSE./size(Y_nanfree,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate VAR with factor structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_OLS = chol(SIGMA_OLS\eye(size(Y_nanfree,2)));
[Y_fact,~] = simvardgp_DP(T,size(Y_nanfree,2),p,A_OLS,PSI_OLS,1);

% Make sure the VAR is stationary
% need to rewrite the VAR(p) as a VAR(1) before checking for stationarity
N = size(Y_nanfree,2);
F = zeros(N,N*p);
for i=1:p
    F(:,1+(i-1)*N:+i*N) = A_OLS(2+(i-1)*N:1+i*N,:);
end
F = [F;[kron(eye(N),eye(p-1)),zeros(N*(p-1),N)]];

lam = eig(F);

for j=1:length(lam)
    if isreal(lam(j))
        check_stat(j)  = abs(lam(j)) >= 1;
    else
        check_stat(j)  = sqrt(real(lam(2))^2 + imag(lam(2))^2) >= 1;
    end
end

if max(check_stat) > 0
    warning('You have generated a nonstationary VAR, please choose different PHI')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate factor loadings and use them to simulate DGP with factor
% structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loadings are generated from uniform (-1,1)
lam_0 = -1 + 2.*rand(K,1);
lam_1 = -1 + 2.*rand(K,size(Y_fact,2));

Y_sim =NaN(T,K);
for t = 1:T
    u = chol(SIGMA)'*randn(K,1);
    Y_sim(t,:) = lam_0 + lam_1 * Y_fact(t,:)' + u;
end




