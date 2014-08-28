function Gamma = envPoisson(Y,X,alpha,gamma,eta)
p = size(X,2);
N = size(X,1);
beta = gamma*eta;
u = size(gamma,2);

%--- get sample statistics ................................................
data_parameters.Y = Y;
data_parameters.X = X;
data_parameters.alpha = alpha;
data_parameters.beta = beta;
data_parameters.eta = eta;
data_parameters.Sx = cov(X)*N/(N-1);
[fooM fooU] = Poisson_cov(Y,X,alpha,beta);
data_parameters.Sxw = fooM;



%--- get handle to objective function and derivative ......................
Fhandle = F(@F4Poisson,data_parameters);
dFhandle = dF(@dF4Poisson,data_parameters);
W0 = gamma;
[fn Gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:u},'quiet',200);