function [dx,dy] = lrt(Y,X,mindim,varargin)
% G is the basis for the sigma_X envelope
% H is the basis for the sigma_Y envelope
% for fixed G, the objective function is the same as the objective function
% of envelope model of Y on G'*X
% for fixed H, the objective function is the same as the objective function
% of envelope model of X on H'*Y

G = eye(size(X,2));
H = eye(size(Y,2));
numiter = 2;
for k=1:numiter
    [HX,H,L,dhat] = mlm_fit(Y,X*G,'lrt',varargin(:));
    [GY,G,L,dhat] = mlm_fit(X,Y*H,'lrt',varargin(:));
    fprintf('k = %d\n',k);
end
dx = max(size(G,2),mindim);
dy = max(size(H,2),mindim);