function [a Gamma e] = Poisson_ini(Y,X,S1,S2,S3)
[N p] = size(X);
data_parameters.Y = Y;
data_parameters.X = X;
data_parameters.Sx = cov(X)*N/(N-1);

data_parameters.alpha = S1.alpha;
data_parameters.eta = S1.eta;
data_parameters.beta = S1.gamma*S1.eta;
[fooM fooU] = Poisson_cov(Y,X,data_parameters.alpha,data_parameters.beta);
data_parameters.Sxw = fooM;
eval1 = F4Poisson(S1.gamma,data_parameters);


data_parameters.alpha = S2.alpha;
data_parameters.eta = S2.eta;
data_parameters.beta = S2.gamma*S2.eta;
[fooM fooU] = Poisson_cov(Y,X,data_parameters.alpha,data_parameters.beta);
data_parameters.Sxw = fooM;
eval2 = F4Poisson(S2.gamma,data_parameters);

data_parameters.alpha = S3.alpha;
data_parameters.eta = S3.eta;
data_parameters.beta = S3.gamma*S3.eta;
[fooM fooU] = Poisson_cov(Y,X,data_parameters.alpha,data_parameters.beta);
data_parameters.Sxw = fooM;
eval3 = F4Poisson(S3.gamma,data_parameters);

Gamma = S1.gamma;
a = S1.alpha;
e = S1.eta;
eval0 = eval1;
if(eval2<eval0)
    Gamma = S2.gamma;
    a = S2.alpha;
    e = S2.eta;
    eval0 = eval2;
end
if(eval3<eval0)
    Gamma = S3.gamma;
    a = S3.alpha;
    e = S3.eta;
    eval0 = eval3;
end
