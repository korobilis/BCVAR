function pvalue=test_CR(benchmark,alternative,H,nested_flag)

T=size(benchmark,1);

differential=benchmark-alternative;

% Pre-weighting
au=differential(2:T,:)\differential(1:T-1,:);
differentialwhite=differential(2:T,:)-differential(1:T-1,:)*au;
T=T-1;

var=NeweyWest(differentialwhite,H);

tdm=sqrt(T)*mean(differentialwhite)/sqrt(var);

% pvalue=1-tcdf(tdm,T);

if nested_flag == 0
    % Models are not nested, 2-sided test. Compute pvalue only if tdm>0,
    % otherwise if tdm<=0, set pvalue to 0
    if tdm > 0
        pvalue=2*tcdf(-abs(tdm),T);
    else
        pvalue = 1;
    end
else
    % Models are nested, 1-sided test
     pvalue=1-tcdf(tdm,T);
end
