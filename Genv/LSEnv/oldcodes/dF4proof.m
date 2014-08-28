function df = dF4proof(W,FParameters)
M1 = FParameters.M1;
M2 = FParameters.M2;
a = M1*W*inv(W'*M1*W);
b = M2*W*inv(W'*M2*W);
df = a + b;
