function Ghat = manifold1D(M,U,d)
% estimating M-envelope contains span(U)
% where M>0 and is symmetric
% dimension of the envelope is d
% based on (M-U) and inv(M)
% if size(U,1)~=size(U,2)
%     U = U*U';
% end
p=size(M,1);
Mnew = M;
Unew = U;
G = zeros(p,d);
G0 = eye(p);
for k=1:d
    gk = first1D(Mnew,Unew);
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Mnew = G0'*M*G0;
    Unew = G0'*U*G0;
end
Ghat = G;