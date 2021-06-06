function [x,mm,ss] = standard(y)
T = size(y,1);
N = size(y,2);

mm = mean(y);
ss = std(y);

my = repmat(mm,T,1);
sy = repmat(ss,T,1);
x = (y-my)./sy;

%x=(y-kron(mean(y),ones(rows(y),1)))./kron(std(y),ones(rows(y),1));
