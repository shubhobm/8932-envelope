function Ghat = envPoisson1D(Y,X,alpha,gamma,eta)
p = size(X,2);
N = size(X,1);
beta = gamma*eta;
u = size(gamma,2);
G = zeros(p,u);
G0 = eye(p);


Ynew = Y;
Xnew = X;


[a_new,G_new,E_new]=Poisson_ini1D(Ynew,Xnew);
for k=1:u
    gk = envPoisson(Ynew,Xnew,a_new,G_new,E_new);
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Xnew = X*G0;
    [a_new,G_new,E_new]=Poisson_ini1D(Ynew,Xnew);
end
Ghat = G;
