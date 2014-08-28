function dim = rrt(Y,X)
% test at level 0.05 the rank of coefficient matrix
% of (reduced rank) regression Y on X
n = size(X,1);
p = size(X,2);
r = size(Y,2);
Bn = inv(X'*X)*X'*Y;
Gn = inv(X'*X./n);
Sigman = (Y-X*Bn)'*(Y-X*Bn)./(n-p-1);
Bstd = inv(sqrtm(Gn))*Bn*inv(sqrtm(Sigman));
minpr = min(p,r);
phi2 = svd(Bstd).^2;
for d = 0:minpr
    Ld = n*sum(phi2((d+1):end));
    df = (r-d)*(p-d);
    pval = 1 - chi2cdf(Ld,df);
    if pval>0.05
        break;
    end
end
dim = d;