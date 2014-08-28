function f = F4cox(W,fpara)
Y = fpara.Y;
X = fpara.X;
delta = fpara.delta;
[N,p] = size(X);
Sxw = fpara.Sxw;
Sxvw = Sxw*fpara.beta;

eta = fpara.eta;
Sx = fpara.Sx;
Xbeta = X*W*inv(W'*Sxw*W)*W'*Sxvw;
vi = zeros(N,1);
for i=1:N
    vi(i) = log((Y>=Y(i))'*exp(Xbeta));
end
Cn = -delta'*Xbeta + delta'*vi;
Mn = N/2*logdet(W'*Sx*W) + N/2*logdet(W'*inv(Sx)*W);
f = Cn + Mn;
