clear all; 
%setpaths;
nrows = 500;
ncols = 20;
nrep = 20;
amax = 15;

% figure 2d
h = 10;
u = 2;
alp = zeros(ncols,2);
alp(1:3,1) = [1 1 1];
alp(1,2)= 1;
alp(ncols-1,2)=1;
alp(ncols,2)=3;
 
angulos = zeros(nrep,amax,5);

for a=1:amax
  disp(strcat('a =',int2str(a)));
  for j=1:nrep
    X=normrnd(0,1,nrows,ncols);
 
%   figure 2d
    yr=.4*a*(X*alp(:,1)).^2 + 3* sin(X*alp(:,2)/4) ;
    y=normrnd(yr,0.2^2);

    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',3);
    angulos(j,a,1)=subspace(W,alp)*180/pi;

    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',5);
    angulos(j,a,2)=subspace(W,alp)*180/pi;
    
    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',10);
    angulos(j,a,3)=subspace(W,alp)*180/pi;

    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',15);
    angulos(j,a,4)=subspace(W,alp)*180/pi;
    
    [WX, W]=fusing(y,X,u);
    angulos(j,a,5)=subspace(W,alp)*180/pi;
  end
end
 
meanang = mean(angulos,1);

plot(squeeze(meanang));
%label
title('Y= X_1/4 + aX_2^2/10 + 3\epsilon/5');
xlabel('a');
ylabel('ANGLE');
legend('h=3','h=5','h=10','h=15','fused','Location','Best');