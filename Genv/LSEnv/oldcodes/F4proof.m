function f = F4proof(W,FParameters)

M1 = FParameters.M1;
M2 = FParameters.M2;
a = log(W' * M1 * W);
b = log(W' * M2 * W);
f = a + b;