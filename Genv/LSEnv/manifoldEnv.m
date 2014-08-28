function Ghat = manifoldEnv(M,U,d,G_ini)
% estimating M-envelope contains span(U)
% where M>0 and is symmetric
% dimension of the envelope is d
% based on (M-U) and inv(M)
if size(U,1)~=size(U,2)
    U = U*U';
end
p = size(M,2);
Diff = M-U+eye(p)*0.0000000000000000000001;
% test if the difference plus a small positive number
% will be positive definite
% equivalently, this is testing if the difference is positive semi-definite
% if not p.s.d, then M=M+U;
[cholA, cholB]=chol(Diff);
if cholB~=0 %Diff is not positive definite
    M = M+U;
end

%--- get sample statistics ................................................
data_parameters.M = M;
data_parameters.U = U;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4manifoldEnv,data_parameters);
dFhandle = dF(@dF4manifoldEnv,data_parameters);
W0 = G_ini;
% p = size(M,2);
% W0 = orth(rand(p,1));
[fn gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',1);
Ghat = gamma;