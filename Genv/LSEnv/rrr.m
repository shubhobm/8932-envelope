function [A B] = rrr( Y, X, k )
% function [A B] = red_rank_reg( Y, X, k )
%
% Solves the system of equations Y= A*B*X + eps where
% A, B are constrained to be Nxk and kxN matrices, i.e., subject
% to rank(AB)=k.
%
% Nick Firoozye
%
 
Sigma=cov([Y,X]);
 
nY = size(Y,2);
nX = size(X,2);
T  = size(X,1);
% = size(Y,1)!!!
 
SigmaYY=Sigma(1:nY,1:nY);
SigmaYX=Sigma(1:nY,nY+1:end);
SigmaXX=Sigma(nY+1:end,nY+1:end);
 
SS = SigmaYX*pinv(SigmaXX)*SigmaYX';
[L,V]=eig(SS);
[U,S,V]=svd(SS);   % eig doesn't work, too low rank!
 
[S ord] = sort(diag(S),'descend');
U = U(:,ord); %reorder columns of V to match newly sorted L
 
A= U(:,1:k);
B = A'*SigmaYX*pinv(SigmaXX);

return;