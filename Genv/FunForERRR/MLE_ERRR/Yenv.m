function [Gamma Eta] = Yenv(X,Y,u,maxitera)
% Envelope MLE
% Y = Gamma*Eta*X + E; dim(Gamma)=r-times-u;
% dim(Eta) = u-times-p; Gamma is semi-orthogonal

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
Gamma = manifoldEnv(Sy,Sxy'*inv(Sx)*Sxy,u,Ghat,maxitera);
Eta = Gamma'*Sxy'*inv(Sx);
