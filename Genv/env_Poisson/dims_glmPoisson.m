function f = dims_glmPoisson(Y,X,alpha)
[n p] = size(X);
N = n;
mux = mean(X);
muy = mean(Y);
Xc = X-mux(ones(N,1),:);
Yc = Y-muy;
M = Xc'*Xc./N;
U = Xc'*Yc*Yc'*Xc./N./(Yc'*Yc);

data_parameters.p = p;
data_parameters.n = n;
data_parameters.N = N;
data_parameters.Sx = M;
data_parameters.Sres = M-U;
data_parameters.Sfit = U;
data_parameters.Y = Y;
data_parameters.X = X;

M = Xc'*Xc./N;
U = Xc'*Yc*Yc'*Xc./N./(Yc'*Yc);

bicval = zeros(p-2,1);
aicval = zeros(p-2,1);
likval = zeros(p-2,1);

for u = 1:min(10,p-2)
    
    [b0 dev0] = glmfit(X,Y,'poisson');
    a0 = b0(1);
    b0 = b0(2:end);

    a_pls = a0;
    b_pls = b0;
    for i=1:3
        [M U] = Poisson_cov(Y,X,a_pls,b_pls);
        G = EnvMU(M,U,u);
        [b2 dev_pls] = glmfit(X*G,Y,'poisson');
        a_pls = b2(1);
        b_pls = G*b2(2:end);
    end
    G_pls = G;

    a_mani = a0;
    b_mani = b0;
    for i=1:3
        [M U] = Poisson_cov(Y,X,a_mani,b_mani);
        G = manifold1D(M,U,u);
        [b2 dev_mani] = glmfit(X*G,Y,'poisson');
        a_mani = b2(1);
        b_mani = G*b2(2:end);
    end
    G_mani = G;
   
    a_lik = a_pls;
    e_lik = G_pls'*b_pls;
    G_lik = G_pls;
    for i=1:3
        G_lik = envPoisson(Y,X,a_lik,G_lik,e_lik);
        [b3 dev_lik] = glmfit(X*G_lik,Y,'poisson');
        a_lik = b3(1);
        e_lik = b3(2:end);
    end
    b_lik = G_lik*e_lik;
    
    LLL = logL_glmPoisson(G_lik,data_parameters,u);
    
    Npar = 1+u+p*(p+1)/2;
    bicval(u) = Npar*log(n)-2*LLL;
    aicval(u) = Npar*2-2*LLL;
    likval(u) = -2*LLL;
end

bicval0=bicval(1);
aicval0 = aicval(1);
bu = 1;
au = 1;
lu = 1;
for u = 2:(p-2)
    if bicval(u)<bicval0
        bu = u;
        bicval0 = bicval(u);
    end
    if aicval(u)<aicval0
        au = u;
        aicval0 = aicval(u);
    end
    chist =  abs(likval(u) - likval(u-1));
    if chist > chi2inv(1-alpha,1)
        lu = u;
    end
end

f = [au bu lu];