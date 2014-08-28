function df = dF4just(W,FParameters)


M = FParameters.M;
U = FParameters.U;
Minv = FParameters.Minv;
a = (M-U)*W*inv(W'*(M-U)*W);
b = Minv*W*inv(W'*Minv*W);
df = a+b;
