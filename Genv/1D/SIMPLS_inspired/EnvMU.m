function Gamma = EnvMU(M,U,d);
p = size(U,1);
W = zeros(p,d);

opts.disp=0;
[W(:,1),ev] = eigs(U,1);
for k=2:d
    Wk = W(:,1:(k-1));
    Ek = M*Wk;
    QEk = eye(p) - Ek*pinv(Ek'*Ek)*Ek';
    [W(:,k),ev] = eigs(QEk*U*QEk,1);
end
Gamma = orth(W);
Gamma0 = null(W');