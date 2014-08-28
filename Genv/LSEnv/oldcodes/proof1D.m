function [fn gamma] = proof1D(M1,M2)
data_parameters.M1 = M1;
data_parameters.M2 = M2;
Fhandle = F(@F4proof,data_parameters);
dFhandle = dF(@dF4proof,data_parameters);
W0 = orth(rand(size(M1,1),1));
[fn gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',500);