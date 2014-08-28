function gamma = first1D(M,U)
% estimating the first 1D direction of the M-envelope contains span(U)
% where M>0 and both M and U are symmetric
p = size(M,2);
% Diff = M-U+eye(p)*0.0000000000000000000001;
% % test if the difference plus a small positive number
% % will be positive definite
% % equivalently, this is testing if the difference is positive semi-definite
% % if not p.s.d, then M=M+U;
% [cholA, cholB]=chol(Diff);
% if cholB~=0 %Diff is not positive definite
%     M = M+U;
% end

%--- get sample statistics ................................................
data_parameters.M = M;
data_parameters.U = U;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4manifold1D,data_parameters);
dFhandle = dF(@dF4manifold1D,data_parameters);
W0 = get_ini1D(M,U);
% p = size(M,2);
% W0 = orth(rand(p,1));
[fn gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',500);