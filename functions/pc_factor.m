function [fac,lam,eigval,ssr] = pc_factor(x,k)

% Compute principal components estimates of factor model
%
% Input:
% x = txn matrix
% k = number of factors
%
% Model
%
% x(it) = lam(i)'f(t) + u(it)
%
% or
%
% x = f*lam' + u
%
%
% Output:
%
% fac = txk matrix of factors
% lam = t*k matrix of factor loadings
% eigval = kx1 vector of eigenvalues of x*x' (or x'x) ...
% ssr = sum of squared u's
%
% Normalization:
%  fac is normalized so that each column has std dev = 1, thus F'F = t*I(k)
%
% Calculation note:
%  Calculations are speeded by using the smaller of x*x' or x'x


t=size(x,1);
n=size(x,2);

if k > 0
    if n < t
        xx=x'*x;
        [ve,va]=eig(xx);
        [aa,ii] = sort(diag(va),'descend');
        va=aa;
        ve=ve(:,ii);
        eigval=va(1:k);
        lam=ve(:,1:k);
        fac=x*lam;
    else
        xx=x*x';
        [ve,va]=eig(xx);
        [aa,ii] = sort(diag(va),'descend');
        va=aa;
        ve=ve(:,ii);
        eigval=va(1:k);
        fac=ve(:,1:k);
        lam=x'*fac;
    end
    
    % Normalize 
    sfac=sqrt(mean(fac.^2));
    fac=fac./repmat(sfac,t,1);
    lam=(x'*fac)/t;  % Note fac'fac = t 
    ssr=sum(va)-sum(eigval);
else
    ssr=sum(sum(x.^2));
    lam=NaN;
    eigval=NaN;
    fac=NaN;
end

