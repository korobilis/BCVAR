function [Q] = grams(A)
[m,~] = size(A);
Q=0*A;
% compute QR using Gram-Schmidt
for j = 1:m
   v = A(j,:);
   for i=1:j-1
        v = v - Q(i,:)*A(j,:)'*Q(i,:);
   end
   Q(j,:) = v/norm(v);
end

