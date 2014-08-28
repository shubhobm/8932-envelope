function f = F4Logistic(W,fpara)
Y = fpara.Y;
X = fpara.X;
N = size(X,1);
Sxw = fpara.Sxw;
Sxvw = Sxw*fpara.beta;

alpha = fpara.alpha;
eta = fpara.eta;
Sx = fpara.Sx;
theta = alpha + X*W*inv(W'*Sxw*W)*W'*Sxvw;
Cn = - Y'*theta + sum(log(1+exp(theta)));
Mn = N/2*logdet(W'*Sx*W) + N/2*logdet(W'*inv(Sx)*W);
f = Cn + Mn;
