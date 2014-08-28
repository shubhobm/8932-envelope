X=fearn(:,1);
Y=fearn(:,2:7);
[GX,G,L,dhat] = mlm_fit(Y,X,'bic');
[betaem, eta, Omega, Omega0, s1, s2]  =  mlm_empars(Y,X,G);
mlm_seratios(Y,X,G);
mlm_fmses(Y,X)
mlm_emses(Y,X,G)