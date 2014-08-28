function f = F4manifold1D(W,FParameters)

M = FParameters.M;
U = FParameters.U;
a = log(W' * (M-U) * W);
b = log(W' * inv(M) * W);
f = a + b;