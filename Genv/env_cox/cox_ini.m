function [Gamma E] = cox_ini(Y,X,delta,S1,S2,S3)
[N p] = size(X);
data_parameters.Y = Y;
data_parameters.X = X;
data_parameters.delta = delta;
data_parameters.Sx = cov(X)*N/(N-1);


data_parameters.eta = S1.eta;
data_parameters.beta = S1.gamma*S1.eta;
[fooM fooU] = cox_cov(Y,X,delta,data_parameters.beta);
data_parameters.Sxw = fooM;
eval1 = F4cox(S1.gamma,data_parameters);


data_parameters.eta = S2.eta;
data_parameters.beta = S2.gamma*S2.eta;
[fooM fooU] = cox_cov(Y,X,delta,data_parameters.beta);
data_parameters.Sxw = fooM;
eval2 = F4cox(S2.gamma,data_parameters);


data_parameters.eta = S3.eta;
data_parameters.beta = S3.gamma*S3.eta;
[fooM fooU] = cox_cov(Y,X,delta,data_parameters.beta);
data_parameters.Sxw = fooM;
eval3 = F4cox(S3.gamma,data_parameters);

Gamma = S1.gamma;
E = S1.eta;
eval0 = eval1;
if(eval2<eval0)
    Gamma = S2.gamma;
    E = S2.eta;
    eval0 = eval2;
end
if(eval3<eval0)
    Gamma = S3.gamma;
    E = S3.eta;
    evl = eval3;
end

