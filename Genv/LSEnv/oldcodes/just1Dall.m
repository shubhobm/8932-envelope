function Ghat = just1Dall(M,U,Minv,d)
p=size(M,1);
Mnew = M;
Unew = U;
Minvnew = Minv;
G = zeros(p,d);
G0 = eye(p);
for k=1:d
    gk = just1D(Mnew,Unew,Minvnew);
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Mnew = G0'*M*G0;
    Unew = G0'*U*G0;
    Minvnew = G0'*inv(M)*G0;
end
Ghat = G;