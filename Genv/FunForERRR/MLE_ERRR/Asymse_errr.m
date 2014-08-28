function [ase_ols ase_rrr ase_env ase_errr] = Asymse_errr(X,Y,d,u,maxitera)
[N p] = size(X);
r = size(Y,2);
datavec = [X Y];
mu = mean(datavec);
datavec = datavec - mu(ones(N,1),:);
Xc = datavec(:,1:p);
Yc = datavec(:,(p+1):(p+r));
Sc = datavec'*datavec./N;
SigmaX = Sc(1:p,1:p);
SigmaY = Sc((p+1):(p+r),(p+1):(p+r));
SigmaXY = Sc(1:p,(p+1):(p+r));
SigmaYX = SigmaXY';
Sigma = SigmaY - SigmaYX*inv(SigmaX)*SigmaXY;


[Aml Bi Ci Bml Cml] = YenvRRR(Xc,Yc,u,d,maxitera);
Gamma = Aml;
eta = Bml;
B = Cml';
A = Gamma*eta;
Gamma0 = null(Gamma');
Omega0 = Gamma0'*Sigma*Gamma0;
Omega = Gamma'*(SigmaY-SigmaYX*B'*B*SigmaXY)*Gamma;
Sigma = Gamma*Omega*Gamma' + Gamma0*Omega0*Gamma0';

%%   OLS
avar_ols = kron(inv(SigmaX),Sigma);

%%   RRR
MA = A*inv(A'*inv(Sigma)*A)*A';
MB = B'*inv(B*SigmaX*B')*B;
avar_rrr = kron(inv(SigmaX),MA) + kron(MB,Sigma) - kron(MB,MA);

%%  ENV
Omega = Gamma'*Sigma*Gamma;
eeta = eta*B;
M1 = kron(eeta*SigmaX*eeta',inv(Omega0))+...
    kron(Omega,inv(Omega0))+...
    kron(inv(Omega),Omega0) - 2*eye(((r-u)*u));
avar_env = kron(inv(SigmaX),Gamma*Omega*Gamma')+...
    kron(eeta',Gamma0)*pinv(M1)*kron(eeta,Gamma0');

%%  ERRR
Omega = Gamma'*(SigmaY-SigmaYX*B'*B*SigmaXY)*Gamma;
Meta = eta*inv(eta'*inv(Omega)*eta)*eta';
Meta = eta*inv(eta'*inv(Omega)*eta)*eta';
avar_errr = kron(MB,Gamma*Omega*Gamma')+...
    kron(inv(SigmaX)-MB,Gamma*Meta*Gamma')+...
    kron(eeta',Gamma0)*pinv(M1)*kron(eeta,Gamma0');

ase_ols = sqrt(diag(avar_ols))./sqrt(N);
ase_env = sqrt(diag(avar_env))./sqrt(N);
ase_rrr = sqrt(diag(avar_rrr))./sqrt(N);
ase_errr = sqrt(diag(avar_errr))./sqrt(N);

