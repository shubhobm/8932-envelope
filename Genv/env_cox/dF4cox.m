function f = dF4cox(W,fpara)
Y = fpara.Y;
X = fpara.X;
delta = fpara.delta;
Sx = fpara.Sx;
[N,p] = size(X);

Sxw = fpara.Sxw;
Sxvw = Sxw*fpara.beta;

Iw = inv(W'*Sxw*W);
eta = Iw*W'*Sxvw;
Xbeta = X*W*eta;

Vn = zeros(N,p);
for i=1:N
	Vn(i,:) = X'*((Y>=Y(i)).*exp(Xbeta))/((Y>=Y(i))'*exp(Xbeta));
end
XC = -X'*delta + Vn'*delta;

M1 = XC*eta';
M2 = Sxvw*XC'*W*inv(W'*Sxw*W);
M3 = -Sxw*W*Iw*W'*XC*eta';
M4 = -Sxw*W*eta*XC'*W*Iw;

Cn = M1+M2+M3+M4;
Mn = N*Sx*W*inv(W'*Sx*W) + N*inv(Sx)*W*inv(W'*inv(Sx)*W);
f = Cn + Mn;