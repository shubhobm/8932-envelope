function [G,H] = AA_env(Yc,Xc,m,l,Sc,Sd)
G = zeros(size(Xc,2),l);
G0 = eye(size(Xc,2));
H = zeros(size(Yc,2),m);
H0 = eye(size(Yc,2));
Xc0 = Xc;
Yc0 = Yc;
for k=1:l
    [GY,gk,L,dhat] = mlm_fit(Xc,Yc0,1); 
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Xc = Xc0*G0;
end
for k=1:m
    [HX,hk,L,dhat] = mlm_fit(Yc,Xc0,1);
    H(:,k) = H0*hk;
    H0 = null(H(:,1:k)');
    Yc = Yc0*H0;
end
G = orth(G);
H = orth(H);
% 
% [HX,H1,L,dhat] = mlm_fit(Yc0,Xc0*G,m,'initval',H);
% [GY,G1,L,dhat] = mlm_fit(Xc0,Yc0*H,l,'initval',G);

% Lik1 = logL(Sc,G1,H1,100);
% Lik = logL(Sc,G,H,100);
% 
% if Lik1<Lik
%     G = G1;
%     H = H1;
% end