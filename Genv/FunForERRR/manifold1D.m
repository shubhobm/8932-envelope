function Ghat = manifold1D(M,U,d)
% estimating M-envelope contains span(U)
% where M>0 and is symmetric
% dimension of the envelope is d
% based on (M-U) and inv(M)
% if size(U,1)~=size(U,2)
%     U = U*U';
% end
if size(U,1)~=size(U,2)
    U = U*U';
end
p = size(M,2);
Diff = M-U+eye(p)*0.0000000000000000000001;
Diff = Diff + Diff';
% test if the difference plus a small positive number
% will be positive definite
% equivalently, this is testing if the difference is positive semi-definite
% if not p.s.d, then M=M+U;
[cholA, cholB]=chol(Diff);
if cholB~=0 %Diff is not positive definite
    M = M+U;
end

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