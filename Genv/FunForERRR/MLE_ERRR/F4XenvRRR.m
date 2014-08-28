function f = F4XenvRRR(W,FParameters)
u = FParameters.u;
d = FParameters.d;
r = FParameters.r;
Sx = FParameters.Sx;
Sres = FParameters.Sres;
aux=inv(sqrtm(W'*Sres*W)) * (W'*Sx*W)*inv(sqrtm(W'*Sres*W));
[vv,aux1]=firsteigs(aux,1); 
f = log(aux1);
a = logdet(W'*Sx*W);
b = logdet(W'*inv(Sx)*W);
aux=inv(sqrtm(W'*Sres*W)) * (W'*Sx*W)*inv(sqrtm(W'*Sres*W));
[vv,aux1]=firsteigs(aux,d); 
f = a + b - sum(log(aux1));

% a = logdet(W'*Sres*W);
% b = logdet(W'*inv(Sx)*W);
% f = a + b;
% 
% Sfit = Sx - Sres;
% aux=inv(sqrtm(W'*Sres*W))*(W'*Sfit*W)*inv(sqrtm(W'*Sres*W));
% [vv,aux1]=firsteigs(aux,u); 
% if d<min(u,r)
%     for i = (d+1):min(u,r)
%         f = f + log(1+aux1(i));
%     end
% end
