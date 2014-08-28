function f = logL_glmPoisson(W,FParameters,u)
p = FParameters.p;
n = FParameters.n;
N = n;
Sx = FParameters.Sx;
Sres = FParameters.Sres;
Sfit = FParameters.Sfit;
Y = FParameters.Y;
X = FParameters.X;

[bw dev_lik] = glmfit(X*W,Y,'poisson');
alik = bw(1);
blik = W*bw(2:end);
theta = alik + X*blik;
Cn = - Y'*theta + sum(exp(theta));
Mn = N/2*logdet(W'*Sx*W) + N/2*logdet(W'*inv(Sx)*W);
f = -Cn - Mn;