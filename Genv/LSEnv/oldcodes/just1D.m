function gamma = just1D(M,U,Minv)
data_parameters.U = U;
data_parameters.M = M;
data_parameters.Minv = Minv;
Fhandle = F(@F4just,data_parameters);
dFhandle = dF(@dF4just,data_parameters);
W0 = get_ini1D(M,U);
[fn gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',500);