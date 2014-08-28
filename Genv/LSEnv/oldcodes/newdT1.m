clear all; 
rand('state',1986)
randn('state',28960)
setpaths;
p = 5; % number of predictors
r = 6; % number of responses
l = 1; % dimension of X-envelope
m = 2; % dimension of Y-envelope
N = 200; % number of observations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Envelopes basis          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueG = orth(rand(p,l));
trueG0 = null(trueG');
trueH = orth(rand(r,m));
trueH0 = null(trueH');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Model parameters         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eta = 10*rand(l,m);
Omega = rand(l,l);
Omega = 15*Omega*Omega';
Omega0 = rand(p-l,p-l);
Omega0 = 5*Omega0*Omega0';
Phi = rand(m,m);
Phi = 5*Phi*Phi';
Phi0 = rand(r-m,r-m);
Phi0 = 10*Phi0*Phi0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Covariance matrices      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaXY = trueG*Omega*Eta*trueH';
sigmaX = trueG*Omega*trueG' + trueG0*Omega0*trueG0';
sigmaY = trueH*(Phi+Eta'*Omega*Eta)*trueH' + trueH0*Phi0*trueH0';
sigmaC = [sigmaX, sigmaXY; sigmaXY', sigmaY];
sigmaD = blkdiag(sigmaX,sigmaY);
sigmaOff = sigmaC - sigmaD;

beta = trueG*Eta*trueH';



%%%%%%%%%%%%%%%%%%%%%%
nsim = 5;  %number of iterations
ang_c = zeros(nsim,2);
ang_xy = zeros(nsim,2);
ang_env = zeros(nsim,2);
ang_pls = zeros(nsim,2);

dev_pls = zeros(nsim,1);
dev_ols = zeros(nsim,1);
dev_c = zeros(nsim,1);
dev_xy = zeros(nsim,1);
dev_env = zeros(nsim,1);
dev_apls = zeros(nsim,1);
Lik_fm = zeros(nsim,1);
Lik_env = zeros(nsim,1);
dx_aic = zeros(nsim,1);
dy_aic = zeros(nsim,1);
dx_bic = zeros(nsim,1);
dy_bic = zeros(nsim,1);
dx_lrt = zeros(nsim,1);
dy_lrt = zeros(nsim,1);
dxy = zeros(nsim,1);

for isim=1:nsim
fprintf('isim = %d\n',isim);
datavec = mvnrnd(zeros(p+r,1),sigmaC,N);
mu = mean(datavec);
datavec = datavec - mu(ones(N,1),:);
Xc = datavec(:,1:p);
Yc = datavec(:,(p+1):(p+r));


datavecN = mvnrnd(zeros(p+r,1),sigmaC,N);
muN = mean(datavecN);
datavecN = datavecN - muN(ones(N,1),:);
XcN = datavecN(:,1:p);
YcN = datavecN(:,(p+1):(p+r));

Sc = datavec'*datavec./N;
Sx = Sc(1:p,1:p);
Sy = Sc((p+1):(p+r),(p+1):(p+r));
Sxy = Sc(1:p,(p+1):(p+r));
Sd = blkdiag(Sx,Sy);
Soff = Sc - Sd;

ScN = datavecN'*datavecN./N;
SxN = ScN(1:p,1:p);
SxyN = ScN(1:p,(p+1):(p+r));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Cook's envelope algorithm       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% OLS %%%%%%%%%%%%%%%%%%%%%%%
betahat = inv(Sx)*Sxy;
resi = XcN*(betahat-beta);
dev_ols(isim) =  sqrt(trace(resi'*resi*resi'*resi));

%%%%%%%%%%%%%%%%%% PLS %%%%%%%%%%%%%%%%%%%%%%%
[XL,YL,XS,YS, betahat] = plsregress(Xc,Yc,max(l,m));
betahat = betahat(2:end,:);
resi = XcN*(betahat-beta);
dev_pls(isim) = sqrt(trace(resi'*resi*resi'*resi));

% %%%%%%%%%% X, Y separately %%%%%%%%%%%%%%%%%%%
% Ghat = EnvMU(Sx,Sxy,l);
% Ghat = orth(Ghat);
% %Hhat = EnvMU(Sy,Sxy',m);
% Hhat = EnvMU(sigmaY,sigmaXY',m);
% Hhat = orth(Hhat);
% ang_xy(isim,1)=subspace(Ghat,trueG)*180/pi;
% ang_xy(isim,2)=subspace(Hhat,trueH)*180/pi;
% betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
% resi = XcN*(betahat-beta);
% dev_xy(isim) =  sqrt(trace(resi'*resi*resi'*resi));
% 
% 



%%%%%%%%%% simultaneously %%%%%%%%%%%%%%%%%%%%
Khat = EnvMU(Sd,Soff,l+m);
%Khat = EnvMU(sigmaD,sigmaOff,l+m);
Ghat = orth(Khat(1:p,:));
Hhat = orth(Khat((p+1):(p+r),:));
ang_c(isim,1)=subspace(Ghat,trueG)*180/pi;
ang_c(isim,2)=subspace(Hhat,trueH)*180/pi;
betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
resi = XcN*(betahat-beta);
dev_c(isim) =  sqrt(trace(resi'*resi*resi'*resi));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Alternating PLS algorithm          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ghat,Hhat] = alter_pls(Yc,Xc,m,l);
ang_pls(isim,1)=subspace(Ghat,trueG)*180/pi;
ang_pls(isim,2)=subspace(Hhat,trueH)*180/pi;
betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
resi = XcN*(betahat-beta);
dev_apls(isim) =  sqrt(trace(resi'*resi*resi'*resi));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Alternating envelope algorithm     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ghat,Hhat] = alter_env(Yc,Xc,m,l,Sc,Sd,Ghat,Hhat);
ang_env(isim,1)=subspace(Ghat,trueG)*180/pi;
ang_env(isim,2)=subspace(Hhat,trueH)*180/pi;
betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
resi = XcN*(betahat-beta);
dev_env(isim) =  sqrt(trace(resi'*resi*resi'*resi));
Lik_fm(isim) = logL(Sc,eye(p),eye(r),N);
Lik_env(isim) = logL(Sc,Ghat,Hhat,N);
% 
 dxy(isim) = rrt(Yc,Xc);
 [dx_lrt(isim), dy_lrt(isim)]=lrt(Yc,Xc,dxy(isim));
 [dx_aic(isim), dy_aic(isim)]=aict(Yc,Xc,dxy(isim),Sc,Sd,N);
 [dx_bic(isim), dy_bic(isim)]=bict(Yc,Xc,dxy(isim),Sc,Sd,N);

end



for dims = 1:5
    fprintf('Dimension=%d\n',dims)
    fprintf('RRT %d\n',sum(dxy==dims))
    fprintf('LRT dx=%d\t dy=%d\n',sum(dx_lrt==dims),sum(dy_lrt==dims))
    fprintf('AIC dx=%d\t dy=%d\n',sum(dx_aic==dims),sum(dy_aic==dims))
    fprintf('BIC dx=%d\t dy=%d\n',sum(dx_bic==dims),sum(dy_bic==dims))
end


% chist = 2*(Lik_env - Lik_fm);

[ang_c, ang_xy, ang_pls,  ang_env]
[dev_ols, dev_pls, dev_c, dev_xy, dev_apls, dev_env]./N

% 
% logdet(Hhat'*sigmaY*Hhat)+logdet(Hhat'*inv(sigmaY)*Hhat)
% 
% logdet(trueH'*sigmaY*trueH)+logdet(trueH'*inv(sigmaY)*trueH)
% 
