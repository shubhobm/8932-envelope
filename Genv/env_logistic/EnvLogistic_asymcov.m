function SEab = EnvLogistic_cov(Y,X,a,Gamma,eta)
N = size(X,1);
p = size(X,2);
u = size(Gamma,2);
b = Gamma*eta;
wts = 1./(2+exp(-a-X*b)+exp(a+X*b));
Sx = cov(X);
Gamma0 = null(Gamma');
Omega = Gamma'*Sx*Gamma;
Omega0 = Gamma0'*Sx*Gamma0;
X = [ones(N,1) X];
V = X'*diag(wts)*X./N;
V = inv(V);
V = V(2:end,2:end);
M = kron(eta,Gamma0')*inv(V)*kron(eta',Gamma0)+kron(Omega,inv(Omega0))+...
    kron(inv(Omega0),Omega)-2*eye(u*(p-u));
V = Gamma*Gamma'*V*Gamma*Gamma'+kron(eta',Gamma0)*inv(M)*kron(eta,Gamma0');
SEab = sqrt(diag(V))./sqrt(N);

% wts = wts./mean(wts);
% Exw = wts'*X/N;
% Sxw = (X-Exw(ones(N,1),:))'*diag(wts)*(X-Exw(ones(N,1),:))./N;
% theta = a + X*b;
% Ys = theta;
% for i=1:N
%     if wts(i)~=0
%         Ys(i) = Ys(i) + (Y(i)-1./(1+exp(-theta(i))))./wts(i);
%     end
% end
% 
% Eyw = wts'*Ys/N;
% Sxyw = (X-Exw(ones(N,1),:))'*diag(wts)*(Ys-Eyw(ones(N,1),:))./N;
% Syw = (Ys-Eyw(ones(N,1),:))'*diag(wts)*(Ys-Eyw(ones(N,1),:))./N;
% M = Sxw;
% U = Sxyw*Sxyw'./Syw;