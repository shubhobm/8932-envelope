function [Ai Bi Aml Bml] = RRR(X,Y,d)
% Y = A*B'*X + E; rank(A)=rank(B)=d;
% A'*A= I_r; A*B'=A*A'*Beta_ols=P_A*Beta_ols;
% two types of solutions: Ai*Bi' are minimizing ||Y-AB'X||
% and Aml*Bml' are the MLE which minimize ||Y-AB'X||_{SigmaYonX}

r = size(Y,2);
p = size(X,2);
n  = size(X,1);
datavec = [Y X];
mu = mean(datavec);
datavec = datavec - mu(ones(n,1),:);
Sigma = datavec'*datavec./n;


SigmaY=Sigma(1:r,1:r);
SigmaYX=Sigma(1:r,r+1:end);
SigmaX=Sigma(r+1:end,r+1:end);

SigmaYonX = SigmaYX*inv(SigmaX)*SigmaYX';
SigmaXonY = SigmaYX'*inv(SigmaY)*SigmaYX;
[L,V]=eig(SigmaYonX);
[U,S,V]=svd(SigmaYonX);   % eig doesn't work, too low rank!
[S ord] = sort(diag(S),'descend');
U = U(:,ord); %reorder columns of V to match newly sorted L
A = U(:,1:d);
B = A'*SigmaYX*inv(SigmaX);
Ai = A;
Bi = B';

% [L,V]=eig(inv(sqrtm(SigmaY))*SigmaYonX*inv(sqrtm(SigmaY)));
% [U,S,V]=svd(inv(sqrtm(SigmaY))*SigmaYonX*inv(sqrtm(SigmaY)));   % eig doesn't work, too low rank!
% [S ord] = sort(diag(S),'descend');
% U = U(:,ord); %reorder columns of V to match newly sorted L
% A = orth(inv(sqrtm(SigmaY))*U(:,1:d));
% B = A'*SigmaYX*inv(SigmaX);
% Aml = A;
% Bml = B';

[L,V]=eig(inv(sqrtm(SigmaX))*SigmaXonY*inv(sqrtm(SigmaX)));
[U,S,V]=svd(inv(sqrtm(SigmaX))*SigmaXonY*inv(sqrtm(SigmaX)));   % eig doesn't work, too low rank!
[S ord] = sort(diag(S),'descend');
U = U(:,ord); %reorder columns of V to match newly sorted L
B = (inv(sqrtm(SigmaX))*U(:,1:d));
A = SigmaYX*B;
Aml = A;
Bml = B;

% [aa,bb] = canoncorr(Y,X);
% Aml = orth(aa(:,1:d));  %r-times-d
% Bml = Aml'*SigmaYX*pinv(SigmaX);
% Bml = Bml';
return;