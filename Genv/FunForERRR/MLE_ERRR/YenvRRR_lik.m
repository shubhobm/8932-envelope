function logLik = YenvRRR_lik(X,Y,u,d,maxitera)
% Envelope MLE
% Y = A*B*C'*X + E; dim(A)=r-times-u;
% dim(B) = u-times-d; dim(C') = d-times-p;
% A and B are semi-orthogonal, so is A*B
[N p] = size(X);
r = size(Y,2);
datavec = [X Y];
mu = mean(datavec);
datavec = datavec - mu(ones(N,1),:);
Xc = datavec(:,1:p);
Yc = datavec(:,(p+1):(p+r));
Sc = datavec'*datavec./N;
Sx = Sc(1:p,1:p);
Sy = Sc((p+1):(p+r),(p+1):(p+r));
Sxy = Sc(1:p,(p+1):(p+r));
Ghat = manifold1D(Sy,Sxy'*inv(Sx)*Sxy,u);
A = XenvRRR(Yc,Xc,u,d,Ghat,maxitera);
W = orth(A);
a = logdet(W'*Sy*W);
b = logdet(W'*inv(Sy)*W);
Sres = Sy - Sxy'*inv(Sx)*Sxy;
aux=inv(sqrtm(W'*Sres*W)) * (W'*Sy*W)*inv(sqrtm(W'*Sres*W));
[vv,aux1]=firsteigs(aux,d); 
f = a + b - sum(log(aux1));
logLik = -N /2*f;