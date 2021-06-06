function out = mvnpdf_DP(y,mean_y,cov_y)

D   = -log(det(cov_y))/2;
Q   = - ((mean_y-y)*inv(cov_y)*(mean_y-y)')/2;
out = exp(-(numel(y)/2)*log(2*pi) + D + Q);