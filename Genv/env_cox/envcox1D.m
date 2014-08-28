function Ghat = envcox1D(Y,X,delta,gamma,eta)
p = size(X,2);
N = size(X,1);
beta = gamma*eta;
u = size(gamma,2);
G = zeros(p,u);
G0 = eye(p);

Ynew = Y;
Xnew = X;

[G_new,E_new]=cox_ini1D(Ynew,Xnew,delta);
for k=1:u
    gk = envcox(Ynew,Xnew,delta,G_new,E_new);
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Xnew = X*G0;
    [G_new,E_new]=cox_ini1D(Ynew,Xnew,delta);
end
Ghat = G;
