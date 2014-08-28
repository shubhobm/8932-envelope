function df = dF4manifoldEnv(W,FParameters)


M = FParameters.M;
U = FParameters.U;
a = (M-U)*W*inv(W'*(M-U)*W);
b = inv(M)*W*inv(W'*inv(M)*W);
df = 2*(a+b);
