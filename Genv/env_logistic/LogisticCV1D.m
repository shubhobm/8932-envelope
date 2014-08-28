function f = LogisticCV1D(Y,X,u,nfold)
[N p] = size(X);
Ntest = floor(N/nfold);
X0 = X;
Y0 = Y;
idx = randperm(N);
Err0 = 0;
Err_pls = 0;
Err_mani = 0;
Err_lik1D = 0;
for i=1:nfold
    testid = [1:Ntest] + (i-1)*Ntest;
    testid = idx(testid);
    X(testid,:) = [];
    Y(testid) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % standard glm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [b0 dev0] = glmfit(X,Y,'binomial','link','logit','constant','on');
    a0 = b0(1);
    b0 = b0(2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMPLS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_pls = a0;
    b_pls = b0;
    for i=1:3
        [M U] = Logistic_cov(Y,X,a_pls,b_pls);
        G = EnvMU(M,U,u);
        [b2 dev_pls] = glmfit(X*G,Y,'binomial','link','logit','constant','on');
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
        [M U] = Logistic_cov(Y,X,a_mani,b_mani);
        G = manifold1D(M,U,u);
        [b2 dev_mani] = glmfit(X*G,Y,'binomial','link','logit','constant','on');
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
    for i=1:3
        G_lik1D = envLogistic1D(Y,X,a_lik1D,G_lik1D,e_lik1D);
        [b3 dev2] = glmfit(X*G_lik1D,Y,'binomial','link','logit','constant','on');
        a_lik1D = b3(1);
        e_lik1D = b3(2:end);
    end
    b_lik1D = G_lik1D*e_lik1D;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xnew = X0(testid,:);
    Ynew = Y0(testid);
    X = X0;
    Y = Y0;
    p0 = exp(a0 + Xnew*b0)./(1+exp(a0 + Xnew*b0));
    p_pls = exp(a_pls + Xnew*b_pls)./(1+exp(a_pls + Xnew*b_pls));
    p_mani = exp(a_mani + Xnew*b_mani)./(1+exp(a_mani + Xnew*b_mani));
    p_lik1D = exp(a_lik1D + Xnew*b_lik1D)./(1+exp(a_lik1D + Xnew*b_lik1D));
    Err0 = Err0 + sum((Ynew - round(p0)).^2);
    Err_pls = Err_pls + sum((Ynew - round(p_pls)).^2);
    Err_mani = Err_mani + sum((Ynew - round(p_mani)).^2);
    Err_lik1D = Err_lik1D + sum((Ynew - round(p_lik1D)).^2);
end

f = [Err0 Err_pls Err_mani Err_lik1D];
