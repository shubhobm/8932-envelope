function f = LogisticCV(Y0,X0,u,nsim,nsample)
[N p] = size(X0);
b0all = zeros(nsim,p);
b_plsall = zeros(nsim,p);
b_maniall = zeros(nsim,p);
b_likall = zeros(nsim,p);

for isim=1:nsim
    idx = randsample(N,nsample,1);
    X = X0(idx,:);
    Y = Y0(idx,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % standard glm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [b0 dev0] = glmfit(X,Y,'binomial','link','logit','constant','on');
    a0 = b0(1);
    b0 = b0(2:end);
    b0all(isim,:) = b0;

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
    b_plsall(isim,:) = b_pls;
    
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
    b_maniall(isim,:) = (G_mani*b2(2:end));


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
    b_likall(isim,:) = b_lik';

end

f = sqrt([diag(cov(b0all)) diag(cov(b_plsall)) diag(cov(b_maniall)) diag(cov(b_likall))])./sqrt(N);
