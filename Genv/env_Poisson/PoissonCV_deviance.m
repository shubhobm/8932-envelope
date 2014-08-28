function f = PoissonCV_deviance(Y,X,u,nfold)
[N p] = size(X);
Ntest = floor(N/nfold);
X0 = X;
Y0 = Y;
idx = randperm(N);
Err0 = 0;
Err_pls = 0;
Err_mani = 0;
Err_lik = 0;
Err_lik1D = 0;
for i=1:nfold
    testid = [1:Ntest] + (i-1)*Ntest;
    testid = idx(testid);
    X(testid,:) = [];
    Y(testid) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % standard glm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [b0 dev0] = glmfit(X,Y,'poisson');
    a0 = b0(1);
    b0 = b0(2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMPLS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D manifold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MLE 1D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_lik1D = a_mani;
    e_lik1D = G_mani'*b_mani;
    G_lik1D = G_mani;
%     for i=1:3
%         G_lik1D = envPoisson1D(Y,X,a_lik1D,G_lik1D,e_lik1D);
%         [b3 dev2] = glmfit(X*G_lik1D,Y,'poisson');
%         a_lik1D = b3(1);
%         e_lik1D = b3(2:end);
%     end
    b_lik1D = G_lik1D*e_lik1D;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     s1.alpha = a_pls;
%     s1.eta = G_pls'*b_pls;
%     s1.gamma = G_pls;
%     s2.alpha = a_mani;
%     s2.eta = G_mani'*b_mani;
%     s2.gamma = G_mani;
%     s3.alpha = a_lik1D;
%     s3.eta = G_lik1D'*b_lik1D;
%     s3.gamma = G_lik1D;
%     [a_lik G_lik e_lik] = Poisson_ini(Y,X,s1,s2,s3);
    
    a_lik = a_mani;
    e_lik = G_mani'*b_mani;
    G_lik = G_mani;
    for i=1:1
        G_lik = envPoisson(Y,X,a_lik,G_lik,e_lik);
        [b3 dev_lik] = glmfit(X*G_lik,Y,'poisson');
        a_lik = b3(1);
        e_lik = b3(2:end);
    end
    b_lik = G_lik*e_lik;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xnew = X0(testid,:);
    Ynew = Y0(testid);
    Xnew(isinf(log(Ynew)),:)=[];
    Ynew(isinf(log(Ynew))) = [];
    X = X0;
    Y = Y0;
    p0 = exp(a0 + Xnew*b0);
    p_pls = exp(a_pls + Xnew*b_pls);
    p_mani = exp(a_mani + Xnew*b_mani);
    p_lik = exp(a_lik + Xnew*b_lik);
    p_lik1D = exp(a_lik1D + Xnew*b_lik1D);
    Err0 = Err0 + sum((log(Ynew) - log(p0)).*p0);
    Err_pls = Err_pls + sum((log(Ynew) - log(p_pls)).*p_pls);
    Err_mani = Err_mani + sum((log(Ynew) - log(p_mani)).*p_mani);
    Err_lik = Err_lik + sum((log(Ynew) - log(p_lik)).*p_lik);
    Err_lik1D = Err_lik1D + sum((log(Ynew) - log(p_lik1D)).*p_lik1D);
end

f = [Err0 Err_pls Err_mani Err_lik Err_lik1D];
