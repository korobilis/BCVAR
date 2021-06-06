function [Zc,K] = create_RHS_cov(Y,M,p,t)

K = M*(M-1)/2;

Zc = zeros((t-p)*M,K);
Stemp = zeros(M,K);
for i = 1:t-p           
    ic = 1;
    for j = 2:M
        Stemp(j,((j-1)+(j-3)*(j-2)/2):ic) = - Y(i+p,1:j-1);
        ic = ic + j;
    end
    Zc((i-1)*M+1:i*M,:) =  Stemp;
end
