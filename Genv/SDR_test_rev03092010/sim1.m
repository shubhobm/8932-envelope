%%%This is the script to reproduce figure 2c from the paper by
%%% R. D. Cook and L. Forzani: "Likelihood-based Sufficient Dimension
%%% Reduction". To appear in JASA
%
% BRIEF DESCRIPTION
% The script compares LAD and F2M methods SIR, SAVE and DR by computing the angle 
% between them and the known central subspace. The regression model for the response 
% is Y = X1/(10*a) + a*X1^2/100 + 3*err/5. Figure shows the average angle for different 
% values of parameter 'a' in the regression model. See the paper for details.
% =========================================================================


clear all; 
%setpaths;
nrows = 500;
ncols = 8;
nrep = 100;
amax = 10;

% figure 2c
 h=5;
 u=1;
 alp = zeros(ncols,1);
 alp(1,1) = 1;

angulos = zeros(nrep,amax,5);

for a=1:amax
  disp(strcat('a =',int2str(a)));
  for j=1:nrep
    X=normrnd(0,1,nrows,ncols);
%   figure 2c
    yr=.6*(X*alp(:,1))/(2*a)+ a*(X*alp(:,1)).^2/(100) ;
    y=normrnd(yr,.6^2);

    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',3);
    angulos(j,a,1)=subspace(W,alp)*180/pi;

    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',5);
    angulos(j,a,2)=subspace(W,alp)*180/pi;
    
    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',10);
    angulos(j,a,3)=subspace(W,alp)*180/pi;

    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',15);
    angulos(j,a,4)=subspace(W,alp)*180/pi;
    
    [WX, W]=fusing(y,X);
    angulos(j,a,5)=subspace(W,alp)*180/pi;
    
  end
end

%meanang = zeros(amax,4);
meanang = mean(angulos,1);

plot(squeeze(meanang));
%label
title('Y=X_1/(10a) + aX_1^2/100+3\epsilon/5');
xlabel('a');
ylabel('ANGLE');
ylim([0 80]);
legend('h=3','h=5','h=10','h=15','fused','Location','Best');