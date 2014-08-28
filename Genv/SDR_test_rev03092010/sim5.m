clear all;
global N;
global p;
p=15;
N=300;
nsim = 100;      


%for(i=1:p)
%    for(j=1:p)
%        sigma(i,j)=0.5^abs(i-j);
%    end
%end
beta=normrnd(0,1,[5,1]);
beta=beta*sqrtm(inv(beta'*beta));
Gamma=normrnd(0,1,[p,2]);
Gamma=Gamma*sqrtm(inv(Gamma'*Gamma));
CS=Gamma;
Pcs=CS*CS'; %projection matrix to the central subspace

hK = 15;    %numers of slices used for lad h=3:hK and kept fusing h=3:15
angles = zeros(hK,nsim);
fangles = zeros(hK,nsim);

for (isim=1:nsim)
    Y=unifrnd(0,5,[N,1]);
    Y=reshape(Y,N,1);
    fY=[(Y<1),(Y>=1&Y<2),(Y>=2&Y<3),(Y>=3&Y<4),(Y>=4&Y<5)];
    gY=[(Y<5/3),(Y>=5/3&Y<10/3),(Y>=10/3&Y<5)];
    X = zeros(N,p);
    for(i=1:N)
        sigmay = fY(i,:)*[1 2 3 4 5]';
        eps = mvnrnd(zeros(p+1,1),eye(p+1),1);
        X(i,:) = eps(1:p)' + sigmay*Gamma(:,2)*eps(p+1);
    end
    X = X + fY*beta*Gamma(:,1)';
    
    for(h=3:hK)
        [WX,W]=ldr(Y,X,'lad','cont',2,'nslices',h);
        angles(h,isim) = subspace(W,CS)*180/pi;
        [WX,W] = fusing(Y,X,2,h);
        fangles(h,isim) = subspace(W,CS)*180/pi;
    end
    disp(isim);
end







>> mean(angles')

ans =

  Columns 1 through 7

         0         0   59.2925   65.7239   55.4005   61.7428   60.2516

  Columns 8 through 14

   62.2007   64.7647   63.1271   65.7153   64.8076   65.8821   66.1138

  Column 15

   68.9096

>> mean(fangles')

ans =

  Columns 1 through 7

         0         0   59.2925   58.0774   53.9885   53.5321   53.0836

  Columns 8 through 14

   54.1813   53.6062   54.2780   53.8284   53.4730   54.1457   52.6885

  Column 15

   54.8225

>> std(angles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    1.4931    1.4073    1.5553    1.5001    1.4466

  Columns 8 through 14

    1.4424    1.3572    1.4707    1.4436    1.3315    1.3161    1.4348

  Column 15

    1.3146

>> std(fangles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    1.4931    1.5137    1.3993    1.4610    1.3336

  Columns 8 through 14

    1.4673    1.4044    1.5151    1.4942    1.3938    1.4933    1.3773

  Column 15

    1.4472
