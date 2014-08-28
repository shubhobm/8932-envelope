function df = dF4XenvRRR(W,FParameters)
% 
% function df = dF4XenvRRR(Y,Fh)
% [N,P] = size(Y);
% ep = 1e-6;
% for k=1:P
%     for j=1:N
%         Yp = Y; Yp(j,k) = Yp(j,k)+ep;
%         Ym = Y; Ym(j,k) = Ym(j,k)-ep;
%         df(j,k) = (Fh(Yp)-Fh(Ym))/ep/2;
%     end
% end


A = FParameters.Sres;
B = FParameters.Sx;
invAw = inv(W'*A*W);
Bw = W'*B*W;
u = FParameters.u;
d = FParameters.d;
Muv = sqrtm(invAw)*Bw*sqrtm(invAw);
[V,Ss,U] = svd(Muv);
V = inv(sqrtm(W'*A*W))*V; % left eigenvector
U = sqrtm(W'*A*W)*U; % right eigenvector
% V = sqrtm(W'*A*W)*V;
% U = sqrtm(invAw)*U;
V = V(:,1:d);
U = U(:,1:d);
T1 = kron(B*W*invAw,eye(u));
T2 = kron(B*W,invAw);
T3 = -kron(A*W*invAw,invAw*Bw);
T4 = -kron(A*W*invAw*Bw,invAw);
df = zeros(size(W));
aux = inv(sqrtm(W'*A*W))*(W'*B*W)*inv(sqrtm(W'*A*W));
[vv,aux1] = firsteigs(Muv,d);
for i=1:d
    lambda_i = aux1(i);
    VecUVt = vec(U(:,i)*V(:,i)');
    VecVUt = vec(V(:,i)*U(:,i)');
    Deriv_i = (T1+T3)*VecVUt + (T2+T4)*VecUVt;
    Deriv_i = reshape(Deriv_i,size(W'))';
    df = df + Deriv_i./lambda_i;
end
a = B*W*inv(W'*B*W);
b = inv(B)*W*inv(W'*inv(B)*W);
df = -df + 2*a + 2*b;




% W = orth(W);
% A = FParameters.Sres;
% B = FParameters.Sx;
% B = B - A;
% invAw = inv(W'*A*W);
% Bw = W'*B*W;
% u = FParameters.u;
% d = FParameters.d;
% r = FParameters.r;
% Muv = sqrtm(invAw)*Bw*sqrtm(invAw);
% [V,Ss,U] = svd(Muv);
% V = inv(sqrtm(W'*A*W))*V; % left eigenvector
% U = sqrtm(W'*A*W)*U; % right eigenvector
% V = V(:,1:u);
% U = U(:,1:u);
% T1 = kron(B*W*invAw,eye(u));
% T2 = kron(B*W,invAw);
% T3 = -kron(A*W*invAw,invAw*Bw);
% T4 = -kron(A*W*invAw*Bw,invAw);
% df = zeros(size(W));
% a = A*W*inv(W'*A*W);
% b = inv(A+B)*W*inv(W'*inv(A+B)*W);
% df = df + a + b;
% 
% if d<min(u,r)
%     aux = inv(sqrtm(W'*A*W))*(W'*B*W)*inv(sqrtm(W'*A*W));
%     [vv,aux1] = firsteigs(aux,u);
%     for i=(d+1):min(u,r)
%         lambda_i = aux1(i);
%         VecUVt = vec(U(:,i)*V(:,i)');
%         VecVUt = vec(V(:,i)*U(:,i)');         
%         Deriv_i = T1*VecUVt + T2*VecVUt + T3*VecUVt + T4*VecUVt;
%         Deriv_i = reshape(Deriv_i,size(W));
%         df = df + Deriv_i./(lambda_i);
%     end
% end


