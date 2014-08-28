function f = dF4Logistic(W,fpara)
Y = fpara.Y;
X = fpara.X;
Sx = fpara.Sx;
N = size(X,1);

Sxw = fpara.Sxw;
Sxvw = Sxw*fpara.beta;

alpha = fpara.alpha;
Iw = inv(W'*Sxw*W);
eta = Iw*W'*Sxvw;
theta = alpha + X*W*eta;
dC = -Y + 1./(1+exp(-theta));
XC = X'*dC;    %p-by-1

M1 = XC*eta';
M2 = Sxvw*XC'*W*inv(W'*Sxw*W);
M3 = -Sxw*W*Iw*W'*XC*eta';
M4 = -Sxw*W*eta*XC'*W*Iw;

Cn = M1+M2+M3+M4;
Mn = N*Sx*W*inv(W'*Sx*W) + N*inv(Sx)*W*inv(W'*inv(Sx)*W);
f = Cn + Mn;