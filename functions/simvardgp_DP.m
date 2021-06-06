function [y,PHI] = simvardgp_DP(T,M,p,PHI,sigma,c_flag)
%--------------------------------------------------------------------------
%   PURPOSE:
%      Get matrix of Y generated from a VAR model
%--------------------------------------------------------------------------
%   INPUTS:
%     T      - Number of observations (rows of Y)
%     M     - Number of series (columns of Y)
%     p      - Number of lags
%     c_flag - 0 for no intercept in the VAR, 1 if intercept should be
%              included
%
%   OUTPUT:
%     y     - [T x N] matrix generated from VAR(p) model
% -------------------------------------------------------------------------

randn('seed',sum(100*clock));
rand('seed',sum(100*clock));

%----------------------GENERATE--------------------------
% Set storage in memory for y
% First p rows are created randomly and are used as 
% starting (initial) values 
y =[rand(p,M) ; zeros(T,M)];

% Now generate Y from VAR (p,PHI,PSI)
for nn = p+1:T+50   
    u = chol(sigma)'*randn(M,1);
    tmp = flipud(y(nn-p:nn-1,:))';
    if c_flag == 0
        y(nn,:) = [tmp(:)']*PHI + u';
    else
        y(nn,:) = [tmp(:)', 1]*PHI + u';
    end
end

y = y(50+1:end,:);

