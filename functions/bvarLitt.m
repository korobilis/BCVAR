function [Ypr,Zpr,Tpr] = bvarLitt(X,p,tau,mu,iRW)


[T,N] = size(X);

% Compute the parameters for the dummies
for jn = 1:N
    [war,Aar,Car]=arfit(X(:,jn),p,p);
    SS0(jn) = sqrt(Car);
end;

MM0 = mean(X);

% Create matrix of regressors and dependent variable
Z = [];

for jp = 1:p
    Z = [Z X(p+1-jp:end-jp,:)];
end

Z = [ones(size(Z(:,1))) Z];

Y = X(p+1:end,:);


% Construct dummy for Litterman prior
Yrw = tau*[zeros(N*(p-1),N); diag(SS0.*iRW)];
Zrw = tau*[zeros(N*p,1) diag(kron(1:p,SS0)) ];

% Construct dummy for the constant
Ycs = 1e-5*[zeros(1,N)];
Zcs = 1e-5*[1 zeros(1,N*p)];

% Construct dummy for the sum of the coefficients
Ylr = mu*diag(MM0.*iRW);
Zlr = mu*[kron(ones(1,p),diag(MM0.*iRW)) zeros(N,1)];

% Construct dummies for prior on covariance matrix of residual;
Ycv = diag(SS0);
Zcv = zeros(N,N*p+1);

% put together all the information
Ypr = [Ycs; Yrw];%; Ylr; Ycv];
Zpr = [Zcs; Zrw];%; Zlr; Zcv];

Tpr = size(Ypr,1);