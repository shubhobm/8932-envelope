function f = proteinCV(Y,X,u,nfold)
[N p] = size(X);
Ntest = floor(N/nfold);
Ntrain = N-Ntest;
X0 = X;
Y0 = Y;
idx = randperm(N);
Err0 = 0;
Err_pls = 0;
Err_mani = 0;
Err_lik1D = 0;
for isim=1:nfold
    testid = [1:Ntest]+isim-1;
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
    % SIMPLS and Mani_1D initialize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Yc = Y - mean(Y);
    mux = mean(X);
    Xc = X - mux(ones(Ntrain,1),:);

    M = Xc'*Xc./Ntrain;
    U = Xc'*Yc*Yc'*Xc./Ntrain./(Yc'*Yc);
    [G_pls,foo] = eigs(U,u);
    e_pls = glmfit(X*G_pls,Y,'binomial','link','logit','constant','on');
    a_pls = e_pls(1);
    b_pls = G_pls*e_pls(2:end);
    G_mani = manifold1D(M,U,u);
    e_mani = glmfit(X*G_mani,Y,'binomial','link','logit','constant','on');
    a_mani = e_mani(1);
    b_mani = G_mani*e_mani(2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MLE 1D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_lik1D = a_mani;
    e_lik1D = G_mani'*b_mani;
    G_lik1D = G_mani;
    for i=1:1
        G_lik1D = envLogistic(Y,X,a_lik1D,G_lik1D,e_lik1D);
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
    p0 = 1./(1+exp(-a0 - Xnew*b0));
    p_pls = 1./(1+exp(-a_pls - Xnew*b_pls));
    p_mani = 1./(1+exp(-a_mani - Xnew*b_mani));
    p_lik1D = 1./(1+exp(-a_lik1D - Xnew*b_lik1D));
    Err0 = Err0 + sum((Ynew - round(p0)).^2);
    Err_pls = Err_pls + sum((Ynew - round(p_pls)).^2);
    Err_mani = Err_mani + sum((Ynew - round(p_mani)).^2);
    Err_lik1D = Err_lik1D + sum((Ynew - round(p_lik1D)).^2);
end

f = [Err0 Err_pls Err_mani Err_lik1D];
