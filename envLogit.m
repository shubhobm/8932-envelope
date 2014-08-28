function outs = envLogit(Y,X,u)
p = size(X,2);
N = size(X,1);
% beta = gamma*eta;
% u = size(gamma,2);

%--- get initial values ........
Yc = Y - mean(Y);
mux = mean(X);
Xc = X - mux(ones(N,1),:);
M = Xc'*Xc./N;
U = Xc'*Yc*Yc'*Xc./N./(Yc'*Yc);
G = manifold1D(M,U,u);
[beta dev_mani] = glmfit(X*G,Y,'binomial','link','logit','constant','on');
alpha = beta(1);
eta = G*beta(2:end);
gamma = G;


%--- get sample statistics ................................................
data_parameters.Y = Y;
data_parameters.X = X;
data_parameters.alpha = alpha;
data_parameters.beta = beta;
data_parameters.eta = eta;
data_parameters.Sx = cov(X)*N/(N-1);
[fooM fooU] = Logistic_cov(Y,X,alpha,beta);
data_parameters.Sxw = fooM;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4Logistic,data_parameters);
dFhandle = dF(@dF4Logistic,data_parameters);
W0 = gamma;
[fn Gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:u},'quiet',125);
outs.Gamma = Gamma;
% outs.Gamma0 = grams(nulbasis(Gamma'));
% outs.beta = Gamma*Gamma'*beta;
% eta = Gamma' * outs.beta;

