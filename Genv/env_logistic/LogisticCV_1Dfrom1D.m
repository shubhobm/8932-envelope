function f = LogisticCV(Y,X,u,nfold)
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
    % MLE 1D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    datavec=[X,Y];
    mu = mean(datavec);
    datavec = datavec - mu(ones(size(datavec,1),1),:);

    Xc = datavec(:,1:p);
    Yc = datavec(:,(p+1):(p+1));

    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1));
    Sy = Sc(p+1,p+1);
    M = Sx;
    U = Sxy*Sxy'./Sy;
    G = manifold1D(M,U,u);
    [b2 dev_mani] = glmfit(X*G,Y,'binomial','link','logit','constant','on');
    a_lik1D = b2(1);
    b_lik1D = G*b2(2:end);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D manifold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_mani = a_lik1D;
    b_mani = b_lik1D;
    for i=1:3
        [M U] = Logistic_cov(Y,X,a_mani,b_mani);
        G = manifold1D(M,U,u);
        [b2 dev_mani] = glmfit(X*G,Y,'binomial','link','logit','constant','on');
        a_mani = b2(1);
        b_mani = G*b2(2:end);
    end
    G_mani = G;



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
    %     [a_lik G_lik e_lik] = Logistic_ini(Y,X,s1,s2,s3);

    a_lik = a_mani;
    e_lik = G_mani'*b_mani;
    G_lik = G_mani;
    for i=1:3
        G_lik = envLogistic(Y,X,a_lik,G_lik,e_lik);
        [b3 dev_lik] = glmfit(X*G_lik,Y,'binomial','link','logit','constant','on');
        a_lik = b3(1);
        e_lik = b3(2:end);
    end
    b_lik = G_lik*e_lik;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xnew = X0(testid,:);
    Ynew = Y0(testid);
    X = X0;
    Y = Y0;
    p0 = 1./(1+exp(-a0 - Xnew*b0));
    p_pls = 1./(1+exp(-a_pls - Xnew*b_pls));
    p_mani = 1./(1+exp(-a_mani - Xnew*b_mani));
    p_lik = 1./(1+exp(-a_lik - Xnew*b_lik));
    p_lik1D = 1./(1+exp(-a_lik1D - Xnew*b_lik1D));
    Err0 = Err0 + sum((Ynew - round(p0)).^2);
    Err_pls = Err_pls + sum((Ynew - round(p_pls)).^2);
    Err_mani = Err_mani + sum((Ynew - round(p_mani)).^2);
    Err_lik = Err_lik + sum((Ynew - round(p_lik)).^2);
    Err_lik1D = Err_lik1D + sum((Ynew - round(p_lik1D)).^2);
end

f = [Err0 Err_pls Err_mani Err_lik Err_lik1D];
