clear all;
rand('state',20130420)
randn('state',20130420)
setpaths;
p = 10; % number of predictors
u = 2; % envelope dimension
N = 50; % number of observations
nsim = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Envelopes basis          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = orth(rand(p,u));
Gamma0 = null(Gamma');
O1 = orth(rand(u,u));
Omega = O1*diag([1 5])*O1';
O2 = orth(rand(p-u,p-u));
%Omega0 = O2*diag([0.02 0.123 0.15*[2:5] 2.765 20])*O2';
Omega0 = 0.2*eye(p-u);
SigmaX = Gamma*Omega*Gamma' + Gamma0*Omega0*Gamma0';
eta = ones(u,1);
beta = Gamma*eta;


vecb = zeros(nsim,4);
vecc = zeros(nsim,4);
vecd = zeros(nsim,4);
for isim = 1:nsim
    X = mvnrnd(zeros(p,1),SigmaX,N);
    Xbeta = X*beta;
    Y =  wblrnd(exp(Xbeta/5),5*ones(N,1));
    delta = floor(2*rand(N,1));
    b0 = coxphfit(X,Y,'censoring',delta);
    dev0 = sqrt((beta-b0)'*(beta-b0));

    b_pls = b0;    
    for i=1:3
        [M U] = cox_cov(Y,X,delta,b_pls);
        G = EnvMU(M,U,u);
        b_pls = coxphfit(X*G,Y,'censoring',delta);
        b_pls = G*b_pls;
    end
    dev_pls = sqrt((beta-b_pls)'*(beta-b_pls));
    G_pls = G;
    
    
    b_mani = b0;    
    for i=1:3
        [M U] = cox_cov(Y,X,delta,b_pls);
        G = manifold1D(M,U,u);
        b_mani = coxphfit(X*G,Y,'censoring',delta);
        b_mani = G*b_mani;
    end
    dev_mani = sqrt((beta-b_mani)'*(beta-b_mani));
    G_mani = G;

    
    G_lik = G_mani;
    e_lik = G_mani'*b_mani;
    for i=1:3
        G_lik = envcox(Y,X,delta,G_lik,e_lik);
        e_lik = coxphfit(X*G_lik,Y,'censoring',delta);
    end
    b_lik = G_lik*e_lik;
    dev_lik = sqrt((beta-b_lik)'*(beta-b_lik));
    
    vecb(isim,:) = [subspace(beta,b0) subspace(beta,b_pls) subspace(beta,b_mani) subspace(beta,b_lik)];
    vecd(isim,:) = [dev0 dev_pls dev_mani dev_lik];
end


vecb = vecb.*180/pi;
mean(vecb)
mean(vecd)
std(vecb)
std(vecd)





clear all;
rand('state',20130420)
randn('state',20130420)
setpaths;
p = 10; % number of predictors
u = 2; % envelope dimension
N = 200; % number of observations
nsim = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Envelopes basis          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = orth(rand(p,u));
Gamma0 = null(Gamma');
O1 = orth(rand(u,u));
Omega = O1*diag([1 5])*O1';
O2 = orth(rand(p-u,p-u));
%Omega0 = O2*diag([0.02 0.123 0.15*[2:5] 2.765 20])*O2';
Omega0 = 0.2*eye(p-u);
SigmaX = Gamma*Omega*Gamma' + Gamma0*Omega0*Gamma0';
eta = ones(u,1);
beta = Gamma*eta;


vecb = zeros(nsim,4);
vecc = zeros(nsim,4);
vecd = zeros(nsim,4);
for isim = 1:nsim
    X = mvnrnd(zeros(p,1),SigmaX,N);
    Xbeta = X*beta;
    Y =  wblrnd(exp(Xbeta/5),5*ones(N,1));
    delta = floor(2*rand(N,1));
    b0 = coxphfit(X,Y,'censoring',delta);
    dev0 = sqrt((beta-b0)'*(beta-b0));

    b_pls = b0;    
    for i=1:3
        [M U] = cox_cov(Y,X,delta,b_pls);
        G = EnvMU(M,U,u);
        b_pls = coxphfit(X*G,Y,'censoring',delta);
        b_pls = G*b_pls;
    end
    dev_pls = sqrt((beta-b_pls)'*(beta-b_pls));
    G_pls = G;
    
    
    b_mani = b0;    
    for i=1:3
        [M U] = cox_cov(Y,X,delta,b_pls);
        G = manifold1D(M,U,u);
        b_mani = coxphfit(X*G,Y,'censoring',delta);
        b_mani = G*b_mani;
    end
    dev_mani = sqrt((beta-b_mani)'*(beta-b_mani));
    G_mani = G;

    
    G_lik = G_mani;
    e_lik = G_mani'*b_mani;
    for i=1:3
        G_lik = envcox(Y,X,delta,G_lik,e_lik);
        e_lik = coxphfit(X*G_lik,Y,'censoring',delta);
    end
    b_lik = G_lik*e_lik;
    dev_lik = sqrt((beta-b_lik)'*(beta-b_lik));
    
    vecb(isim,:) = [subspace(beta,b0) subspace(beta,b_pls) subspace(beta,b_mani) subspace(beta,b_lik)];
    vecd(isim,:) = [dev0 dev_pls dev_mani dev_lik];
end


vecb = vecb.*180/pi;
mean(vecb)
mean(vecd)
std(vecb)
std(vecd)





clear all;
rand('state',20130420)
randn('state',20130420)
setpaths;
p = 10; % number of predictors
u = 2; % envelope dimension
N = 800; % number of observations
nsim = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Envelopes basis          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = orth(rand(p,u));
Gamma0 = null(Gamma');
O1 = orth(rand(u,u));
Omega = O1*diag([1 5])*O1';
O2 = orth(rand(p-u,p-u));
%Omega0 = O2*diag([0.02 0.123 0.15*[2:5] 2.765 20])*O2';
Omega0 = 0.2*eye(p-u);
SigmaX = Gamma*Omega*Gamma' + Gamma0*Omega0*Gamma0';
eta = ones(u,1);
beta = Gamma*eta;


vecb = zeros(nsim,4);
vecc = zeros(nsim,4);
vecd = zeros(nsim,4);
for isim = 1:nsim
    X = mvnrnd(zeros(p,1),SigmaX,N);
    Xbeta = X*beta;
    Y =  wblrnd(exp(Xbeta/5),5*ones(N,1));
    delta = floor(2*rand(N,1));
    b0 = coxphfit(X,Y,'censoring',delta);
    dev0 = sqrt((beta-b0)'*(beta-b0));

    b_pls = b0;    
    for i=1:3
        [M U] = cox_cov(Y,X,delta,b_pls);
        G = EnvMU(M,U,u);
        b_pls = coxphfit(X*G,Y,'censoring',delta);
        b_pls = G*b_pls;
    end
    dev_pls = sqrt((beta-b_pls)'*(beta-b_pls));
    G_pls = G;
    
    
    b_mani = b0;    
    for i=1:3
        [M U] = cox_cov(Y,X,delta,b_pls);
        G = manifold1D(M,U,u);
        b_mani = coxphfit(X*G,Y,'censoring',delta);
        b_mani = G*b_mani;
    end
    dev_mani = sqrt((beta-b_mani)'*(beta-b_mani));
    G_mani = G;

    
    G_lik = G_mani;
    e_lik = G_mani'*b_mani;
    for i=1:3
        G_lik = envcox(Y,X,delta,G_lik,e_lik);
        e_lik = coxphfit(X*G_lik,Y,'censoring',delta);
    end
    b_lik = G_lik*e_lik;
    dev_lik = sqrt((beta-b_lik)'*(beta-b_lik));
    
    vecb(isim,:) = [subspace(beta,b0) subspace(beta,b_pls) subspace(beta,b_mani) subspace(beta,b_lik)];
    vecd(isim,:) = [dev0 dev_pls dev_mani dev_lik];
end


vecb = vecb.*180/pi;
mean(vecb)
mean(vecd)
std(vecb)
std(vecd)


