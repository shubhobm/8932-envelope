function f = F4just(W,FParameters)

M = FParameters.M;
Minv = FParameters.Minv;
U = FParameters.U;
a = log(W' * (M-U) * W);
b = log(W' * Minv * W);
f = a + b;