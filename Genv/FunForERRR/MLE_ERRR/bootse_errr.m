function [se_ols se_rrr se_env se_errr] = bootse_errr(X,Y,d,u,maxitera,nsim)
[N p] = size(X);
r = size(Y,2);
datavec = [X Y];
mu = mean(datavec);
datavec = datavec - mu(ones(N,1),:);
Xc = datavec(:,1:p);
Yc = datavec(:,(p+1):(p+r));
Sc = datavec'*datavec./N;
Sx = Sc(1:p,1:p);
Sy = Sc((p+1):(p+r),(p+1):(p+r));
Sxy = Sc(1:p,(p+1):(p+r));
Syx = Sxy';

Xc0 = Xc;Yc0 = Yc;Sc0 = Sc;Sxy0 = Sxy;Syx0 = Syx;Sx0=Sx;Sy0=Sy;
beta_ols = zeros(nsim,p*r);
beta_env = zeros(nsim,p*r);
beta_rrr = zeros(nsim,p*r);
beta_errr = zeros(nsim,p*r);
%%%%%%%%%%%%%%%%%% OLS %%%%%%%%%%%%%%%%%%%%%%%
betahat0 = Sxy0'*inv(Sx0);
resi = Yc0-Xc0*betahat0';
for isim = 1:nsim
    Yc = Xc0*betahat0' + resi(randsample(1:N,N,true),:);
    datavec = [Xc0,Yc];
    Sc = datavec'*datavec./N;
    Sx = Sc(1:p,1:p);
    Sxy = Sc(1:p,(p+1):(p+r));
    beta_ols(isim,:) = reshape(Sxy'*inv(Sx),1,p*r);
end

%%%%%%%%%%%%%%%%%% Env %%%%%%%%%%%%%%%%%%%%%%%
[Ghat Ehat] = Yenv(Xc0,Yc0,u,maxitera);
betahat0 = Ghat*Ehat;
resi = Yc0-Xc0*betahat0';
for isim = 1:nsim
    Yc = Xc0*betahat0' + resi(randsample(1:N,N,true),:);
    Xc = Xc0;
    [Ghat Ehat] = Yenv(Xc,Yc,u,maxitera);
    beta_env(isim,:) = reshape(Ghat*Ehat,1,p*r);
end

%%%%%%%%%%%%%%%%%% RRR %%%%%%%%%%%%%%%%%%%%%%%
[Ai Bi Aml Bml] = RRR(Xc0,Yc0,d);
betahat0 = Aml*Bml';
resi = Yc0-Xc0*betahat0';
for isim = 1:nsim
    Yc = Xc0*betahat0' + resi(randsample(1:N,N,true),:);
    Xc = Xc0;
    [Ai Bi Aml Bml] = RRR(Xc0,Yc,d);
    beta_rrr(isim,:) = reshape(Aml*Bml',1,p*r);
end


%%%%%%%%%%%%%%%%%% ERRR %%%%%%%%%%%%%%%%%%%%%%%
[Aml Bi Ci Bml Cml] = YenvRRR(Xc0,Yc0,u,d,maxitera);
betahat0 = Aml*Bml*Cml';
resi = Yc0-Xc0*betahat0';
for isim = 1:nsim
    Yc = Xc0*betahat0' + resi(randsample(1:N,N,true),:);
    Xc = Xc0;
    [Aml Bi Ci Bml Cml] = YenvRRR(Xc0,Yc,u,d,maxitera);
    beta_errr(isim,:) = reshape(Aml*Bml*Cml',1,p*r);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Summary and comparison           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
se_ols = sqrt(diag(get_cov(beta_ols)));
se_env = sqrt(diag(get_cov(beta_env)));
se_rrr = sqrt(diag(get_cov(beta_rrr)));
se_errr = sqrt(diag(get_cov(beta_errr)));



