function [PRED,beta] = bvarLitt2(X,p,tau,mu,H,iRW);


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

Z = [Z ones(size(Z(:,1)))];

Y = X(p+1:end,:);


% Construct dummy for Litterman prior
Yrw = tau*[diag(SS0.*iRW);zeros(N*(p-1),N)];
Zrw = tau*[diag(kron(1:p,SS0)) zeros(N*p,1)];

% Construct dummy for the constant
Ycs = 1e-5*[zeros(1,N)];
Zcs = 1e-5*[zeros(1,N*p) 1];

% Construct dummy for the sum of the coefficients
Ylr = mu*diag(MM0.*iRW);
Zlr = mu*[kron(ones(1,p),diag(MM0.*iRW)) zeros(N,1)];

% Construct dummies for prior on covariance matrix of residual;
Ycv = diag(SS0);
Zcv = zeros(N,N*p+1);

% put together all the information
Ypr = [Yrw;Ylr; Ycv; Ycs]; Zpr = [Zrw; Zlr; Zcv; Zcs];

Tstar = size(Ypr,1);

% Compute the posterior
ZZinv = inv(Zpr'*Zpr + Z'*Z);
ZY    = Zpr'*Ypr + Z'*Y;

beta = ZZinv*ZY;

SIGMA = (Y-Z*beta)'*(Y-Z*beta)./(size(Y,1)-p);

% Compute the fit
if nargout>1
    SS = 1/(T-p+1)*(Y-Z*beta)'*(Y-Z*beta);
   fit = 1 - diag(SS)./diag(cov(diff(X)));
end;

% Compute the forecasts
pred = [];
for j = 0:H
    
    X = [X;pred];
    Z = [];
    for jp = 0:p-1
        Z = [Z X(end-jp,:)];
    end
    
    Z = [Z 1];
    pred = Z(end,:)*beta;
    
end;

PRED = X(end-H+1:end,:)';