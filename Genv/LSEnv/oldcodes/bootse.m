function [verOLS,verENV] = bootse(datavec,p,r,l,m,N,nsim,dPLS,dCCA);
mu = mean(datavec);
datavec = datavec - mu(ones(N,1),:);
Xc0 = datavec(:,1:p);
Yc0 = datavec(:,(p+1):(p+r));
Sc0 = datavec'*datavec./N;
Sx0 = Sc0(1:p,1:p);
Sy0 = Sc0((p+1):(p+r),(p+1):(p+r));
Sxy0 = Sc0(1:p,(p+1):(p+r));
Sd0 = blkdiag(Sx0,Sy0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betahat0 = inv(Sx0)*Sxy0;
resi = Yc0-Xc0*betahat0;
for isim = 1:nsim
    Yc = Xc0*betahat0 + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    beta_ols(isim,:) = reshape(inv(Sx)*Sxy,1,p*r);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XL,YL,XS,YS, betahat] = plsregress(Xc0,Yc0,dPLS);
betahat0 = betahat(2:end,:);
resi = Yc0-Xc0*betahat0;
for isim = 1:nsim
    Yc = Xc0*betahat0 + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    [XL,YL,XS,YS, betahat] = plsregress(Xc0,Yc,dPLS);
    betahat1 = betahat(2:end,:);    
    beta_pls(isim,:) = reshape(betahat1,1,p*r);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [aa,bb] = canoncorr(Xc0,Yc0);
    Ghat = orth(aa(:,1:dCCA));
    Hhat = orth(bb(:,1:dCCA));
    betahat0 = Ghat*pinv(Ghat'*Sx0*Ghat)*Ghat'*Sxy0*Hhat*Hhat';
    resi = Yc0-Xc0*betahat0;
for isim = 1:nsim
    Yc = Xc0*betahat0 + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    [aa,bb] = canoncorr(Xc0,Yc);
    Ghat = orth(aa(:,1:dCCA));
    Hhat = orth(bb(:,1:dCCA));
    betahat1 = Ghat*pinv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
    beta_cca(isim,:) = reshape(betahat1,1,p*r);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simultaneous Env
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ghat,Hhat] = AA_env(Yc0,Xc0,m,l,Sc0,Sd0);
Ghat = orth(Ghat);
Hhat = orth(Hhat);
betahat0 = Ghat*inv(Ghat'*Sx0*Ghat)*Ghat'*Sxy0*Hhat*Hhat';
resi = Yc0-Xc0*betahat0;
for isim = 1:nsim
    Yc = Xc0*betahat0 + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    Sy = Sc((p+1):(p+r),(p+1):(p+r));
    Sd = blkdiag(Sx,Sy);
    [Ghat,Hhat] = AA_env(Yc,Xc0,m,l,Sc,Sd);
    Ghat = orth(Ghat);
    Hhat = orth(Hhat);
    betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
    beta_env(isim,:) = reshape(betahat,1,p*r);
    fprintf('%d\n',isim)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   X-Env
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ghat,Hhat] = AA_env(Yc0,Xc0,m,1,Sc0,Sd0);
Ghat = orth(Ghat);
Hhat = orth(Hhat);
Ghat = eye(p);
betahat0 = Ghat*inv(Ghat'*Sx0*Ghat)*Ghat'*Sxy0*Hhat*Hhat';
resi = Yc0-Xc0*betahat0;
for isim = 1:nsim
    Yc = Xc0*betahat0 + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    Sy = Sc((p+1):(p+r),(p+1):(p+r));
    Sd = blkdiag(Sx,Sy);
    [Ghat,Hhat] = AA_env(Yc,Xc0,m,1,Sc,Sd);
    Ghat = orth(Ghat);
    Hhat = orth(Hhat);
    Ghat =  eye(p);
    betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
    beta_Xenv(isim,:) = reshape(betahat,1,p*r);
    fprintf('%d\n',isim)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Y-Env
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ghat,Hhat] = AA_env(Yc0,Xc0,1,l,Sc0,Sd0);
Ghat = orth(Ghat);
Hhat = orth(Hhat);
Hhat = eye(r);
betahat0 = Ghat*inv(Ghat'*Sx0*Ghat)*Ghat'*Sxy0*Hhat*Hhat';
resi = Yc0-Xc0*betahat0;
for isim = 1:nsim
    Yc = Xc0*betahat0 + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    Sy = Sc((p+1):(p+r),(p+1):(p+r));
    Sd = blkdiag(Sx,Sy);
    [Ghat,Hhat] = AA_env(Yc,Xc0,1,l,Sc,Sd);
    Ghat = orth(Ghat);
    Hhat = orth(Hhat);
    Hhat =  eye(r);
    betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
    beta_Yenv(isim,:) = reshape(betahat,1,p*r);
    fprintf('%d\n',isim)
end

verOLS=[[min(sqrt(diag(get_cov(beta_pls)))./sqrt(diag(get_cov(beta_ols)))), max(sqrt(diag(get_cov(beta_pls)))./sqrt(diag(get_cov(beta_ols))))];
[min(sqrt(diag(get_cov(beta_cca)))./sqrt(diag(get_cov(beta_ols)))), max(sqrt(diag(get_cov(beta_cca)))./sqrt(diag(get_cov(beta_ols))))];
[min(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_ols)))), max(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_ols))))];
[min(sqrt(diag(get_cov(beta_Xenv)))./sqrt(diag(get_cov(beta_ols)))), max(sqrt(diag(get_cov(beta_Xenv)))./sqrt(diag(get_cov(beta_ols))))];
[min(sqrt(diag(get_cov(beta_Yenv)))./sqrt(diag(get_cov(beta_ols)))), max(sqrt(diag(get_cov(beta_Yenv)))./sqrt(diag(get_cov(beta_ols))))]]



verENV=[[min(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_ols)))), max(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_ols))))];
[min(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_pls)))), max(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_pls))))];
[min(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_cca)))), max(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_cca))))];
[min(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_Xenv)))), max(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_Xenv))))];
[min(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_Yenv)))), max(sqrt(diag(get_cov(beta_env)))./sqrt(diag(get_cov(beta_Yenv))))]]

