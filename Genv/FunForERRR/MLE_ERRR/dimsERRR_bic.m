function [dbic ubic] = dimsERRR_bic(X,Y,alpha,maxitera)
drrt = RankTest(X,Y,alpha);

n = size(X,1);
p = size(X,2);
r = size(Y,2);
datavec = [X Y];
mu = mean(datavec);
datavec = datavec - mu(ones(n,1),:);
Xc = datavec(:,1:p);
Yc = datavec(:,(p+1):(p+r));

u0 = drrt;
bic1 = zeros(r,3);        % first column d=drrt-1; second column d=drrt; third d=drrt+1;
LogLik0 = YenvRRR_lik(Xc,Yc,u0-1,drrt-1,maxitera);
bic1(u0,1) = (r-u0)*(drrt-1)*log(n) - 2*LogLik0;
LogLik0 = YenvRRR_lik(Xc,Yc,u0,drrt,maxitera);
bic1(u0,2) = (r-u0)*drrt*log(n) - 2*LogLik0;
LogLik0 = YenvRRR_lik(Xc,Yc,u0+1,drrt+1,maxitera);
bic1(u0,3) = (r-u0)*(drrt+1)*log(n) - 2*LogLik0;;    


for u = u0:(r-5)
    LogLik = YenvRRR_lik(Xc,Yc,u,drrt-1,maxitera);
    bic1(u,1) = -((r-u)*(drrt-1)+(p-drrt+1)*(r-drrt+1))*log(n) - 2*LogLik;
end
for u = (u0+1):(r-5)
    LogLik = YenvRRR_lik(Xc,Yc,u,drrt,maxitera);
    bic1(u,2) = -((r-u)*(drrt)+(p-drrt)*(r-drrt))*log(n) - 2*LogLik;
end
for u=(u0+2):(r-5)
    LogLik = YenvRRR_lik(Xc,Yc,u,drrt+1,maxitera);
    bic1(u,3) = -((r-u)*(drrt+1)+(p-drrt-1)*(r-drrt-1))*log(n) - 2*LogLik;
end

[v1 u1] = min(bic1((u0-1):end,1));
[v2 u2] = min(bic1(u0:end,2));
[v3 u3] = min(bic1((u0+1):end,3));
mv = min([v1 v2 v3]);
if v1==mv
    dbic = drrt-1;
    ubic = u1+u0-1;
end
if v2==mv
    dbic = drrt;
    ubic = u2+u0;
end
if v3==mv
    dbic = drrt+1;
    ubic = u3+u0+1;
end

