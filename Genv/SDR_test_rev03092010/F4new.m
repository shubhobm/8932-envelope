function f = F4new(W,FParameters)
% Objective function (minus the likelihood) for the PFC model.
% USAGE: 
% - W is the projection matrix onto the reduced subspace.
% - u is the dimension of the reduced subspace.
% Notice that global FParameters and SEPFCparameters are supposed 
% to be already set.
% ==============================================================
if nargin < 1,
    error('not enough input parameters')
end
u = FParameters.u;
d = FParameters.d;
r = FParameters.r;
A = FParameters.sigmag;
n = sum(FParameters.n);
% p = cols(A);
Afit = FParameters.Afit;
sigres = A-Afit;
% eigtem = eig(W'*sigres*W);
a = logdet(W'*A*W);
% eigtem1 = eig(W'*inv(A)*W);
b =logdet(W'*inv(A)*W);
aux=inv(sqrtm(W'*sigres*W)) * (W'*A*W)*inv(sqrtm(W'*sigres*W));
% [vv,aux1]=firsteigs(aux,u);
% aux11=aux1(d+1:min(u,r));
% f =   n/2*a+n/2*b+n/2*sum(log(1+aux11));
[vv,aux1]=firsteigs(aux,d); 
f =   n/2*a+n/2*b-n/2*sum(log(aux1));