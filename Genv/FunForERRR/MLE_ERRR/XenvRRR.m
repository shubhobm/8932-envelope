function G = XenvRRR(Xc,Yc,u,d,Ghat,maxitera)
p = size(Xc,2);
r = size(Yc,2);
Sc = cov([Xc,Yc]);
%--- get sample statistics ................................................
data_parameters.u = u;
data_parameters.d = d;
data_parameters.r = r;
Sx = Sc(1:p,1:p);
Sxy = Sc(1:p,((p+1):end));
Sy = Sc((p+1):end,(p+1):end);
data_parameters.Sx = Sx;
data_parameters.Sres = Sx-Sxy*inv(Sy)*Sxy';

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4XenvRRR,data_parameters);
%dFhandle = dF(@dF4XenvRRR,Fhandle);
dFhandle = dF(@dF4XenvRRR,data_parameters);
[fn G] = sg_min(Fhandle,dFhandle,Ghat,'prcg','euclidean',{1:u},'quiet',maxitera);
% fini = F4XenvRRR(Ghat,data_parameters);
% ffin = F4XenvRRR(G,data_parameters);
% if fini<ffin
%     G = orth(Ghat);
% end
