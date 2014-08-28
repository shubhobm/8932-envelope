function G = X_env(Xc,Yc,Sc,dim)
p = size(Xc,2);
r = size(Yc,2);
G = zeros(size(Xc,2),dim);
G0 = eye(size(Xc,2));
Xc0 = Xc;
Yc0 = Yc;
for k=1:dim
    [GY,gk,L,dhat] = mlm_fit(Xc,Yc,1); 
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Xc = Xc0*G0;
end
G = orth(G);
% [GY,G1,L,dhat] = mlm_fit(Xc0,Yc0,dim,'initval',G);
% G1 = orth(G1);
% 
% Lik1 = logL(Sc,G1,eye(r),100);
% Lik = logL(Sc,G,eye(r),100);
% if Lik1<Lik
%     G = G1;
% end