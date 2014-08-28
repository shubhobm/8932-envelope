function Gamma = envcox(Y,X,delta,gamma,eta)
[N,p] = size(X);
beta = gamma*eta;
u = size(gamma,2);

%--- get sample statistics ................................................
data_parameters.Y = Y;
data_parameters.X = X;
data_parameters.delta = delta;
data_parameters.beta = beta;
data_parameters.eta = eta;
data_parameters.Sx = cov(X)*N/(N-1);
[fooM fooU] = cox_cov(Y,X,delta,beta);
data_parameters.Sxw = fooM;


%--- get handle to objective function and derivative ......................
Fhandle = F(@F4cox,data_parameters);
dFhandle = dF(@dF4cox,data_parameters);
W0 = gamma;
[fn Gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:u},'quiet',100);