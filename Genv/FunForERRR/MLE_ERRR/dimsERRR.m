function [drrt uaic ubic ulrt] = dimsERRR(X,Y,alpha_d,alpha_u,maxitera)
drrt = RankTest(X,Y,alpha_d);

n = size(X,1);
p = size(X,2);
r = size(Y,2);
datavec = [X Y];
mu = mean(datavec);
datavec = datavec - mu(ones(n,1),:);
Xc = datavec(:,1:p);
Yc = datavec(:,(p+1):(p+r));

u0 = drrt;
LogLik0 = YenvRRR_lik(Xc,Yc,u0,drrt,maxitera);

bic0 = (r-u0)*drrt*log(n) - 2*LogLik0;
bic1 = zeros(r,1);        
bic1(u0) = bic0;    

aic0 = 2*(r-u0)*drrt - 2*LogLik0;
aic1 = zeros(r,1);        
aic1(u0) = aic0;    

Lik0 = 2*LogLik0;
Lik1 = zeros(r,1);
Lik1(u0) = Lik0;

for u = (u0+1):(r-2)
    LogLik = YenvRRR_lik(Xc,Yc,u,drrt,maxitera);
    bic1(u) = -(r-u)*drrt*log(n) - 2*LogLik;
    aic1(u) = -2*(r-u)*drrt - 2*LogLik;
    Lik1(u) = 2*abs(LogLik-LogLik0);
end

[d0 u_bic] = min(bic1(u0:end));
[d0 u_aic] = min(aic1(u0:end));
ubic = u0+u_bic;
uaic = u0+u_aic;

for u = (u0+1):(r-2)
    df = (r-u)*drrt;
    pval = 1 - chi2cdf(Lik1(u),df);
    if pval>alpha_u
        break;
    end
end
ulrt = u;

