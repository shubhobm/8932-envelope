function f = F4manifoldEnv(W,FParameters)

M = FParameters.M;
U = FParameters.U;
a = logdet(W' * (M-U) * W);
b = logdet(W' * inv(M) * W);
f = a + b;