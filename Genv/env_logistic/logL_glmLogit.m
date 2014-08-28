function f = logL_glmLogit(W,FParameters,u)
p = FParameters.p;
n = FParameters.n;
N = n;
Sx = FParameters.Sx;
Sres = FParameters.Sres;
Sfit = FParameters.Sfit;
Y = FParameters.Y;
X = FParameters.X;

[bw dev_lik] = glmfit(X*W,Y,'binomial','link','logit','constant','on');
alik = bw(1);
blik = W*bw(2:end);
theta = alik + X*blik;
Cn = - Y'*theta + sum(log(1+exp(theta)));
Mn = N/2*logdet(W'*Sx*W) + N/2*logdet(W'*inv(Sx)*W);
f = -Cn - Mn;