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
beta=normrnd(0,1,[5,2]);
beta=beta*sqrtm(inv(beta'*beta));
Gamma=normrnd(0,1,[p,2]);
Gamma=Gamma*sqrtm(inv(Gamma'*Gamma));
CS=Gamma;
Pcs=CS*CS'; %projection matrix to the central subspace

hK = 15;    %numers of slices used for lad h=3:hK and kept fusing h=3:15
angles = zeros(hK,nsim);

for (isim=1:nsim)
    Y=unifrnd(0,5,[N,1]);
    Y=reshape(Y,N,1);
    fY=[(Y<1),(Y>=1&Y<2),(Y>=2&Y<3),(Y>=3&Y<4),(Y>=4&Y<5)];
    %gy=[0.1*(Y<1),0.5*(Y>=1&Y<2),0.5*(Y>=2&Y<3),10*(Y>=3&Y<4),100*(Y>=4&Y<5)];
    %gy=sum(gy,2);
    %gy=Y.^2;
    %gy=[gy gy gy gy gy gy];
    X = zeros(N,p);
    for(i=1:N)
        sigmay = Y(i);
        eps = mvnrnd(zeros(p+2,1),eye(p+2),1);
        X(i,:) = eps(1:p)' + sigmay*CS*eps((p+1):(p+2))';
    end
    X = X + 0.5*fY*beta*Gamma';
    
    for(h=3:hK)
        [WX,W]=ldr(Y,X,'lad','cont',2,'nslices',h);
        angles(h,isim) = subspace(W,CS)*180/pi;
    end
    [WX,W] = fusing(Y,X,2,15);
    angles(1,isim)=subspace(W,CS)*180/pi;
    disp(isim);
end











>> sum(angles')/nsim
ans =

  Columns 1 through 9 

   34.5178         0   37.4704   36.0326   35.7767   35.2949   36.0705   36.2510   36.5607

  Columns 10 through 15 

   37.7677   38.2905   39.3119   39.5159   41.3890   42.4737

>> std(angles')/sqrt(nsim)

ans =

  Columns 1 through 9 

    0.4609         0    0.5328    0.4947    0.4638    0.4563    0.4863    0.4881    0.4529

  Columns 10 through 15 

    0.6056    0.5136    0.6044    0.5704    0.7224    0.8411
    

    
    
    
    
    
    
    
    
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
beta=normrnd(0,1,[5,2]);
beta=beta*sqrtm(inv(beta'*beta));
Gamma=normrnd(0,1,[p,2]);
Gamma=Gamma*sqrtm(inv(Gamma'*Gamma));
CS=Gamma;
Pcs=CS*CS'; %projection matrix to the central subspace

hK = 15;    %numers of slices used for lad h=3:hK and kept fusing h=3:15
angles = zeros(hK,nsim);

for (isim=1:nsim)
    disp(isim);
    Y=unifrnd(0,5,[N,1]);
    Y=reshape(Y,N,1);
    fY=[(Y<1),(Y>=1&Y<2),(Y>=2&Y<3),(Y>=3&Y<4),(Y>=4&Y<5)];
    %gy=[0.1*(Y<1),0.5*(Y>=1&Y<2),0.5*(Y>=2&Y<3),10*(Y>=3&Y<4),100*(Y>=4&Y<5)];
    %gy=sum(gy,2);
    %gy=Y.^2;
    %gy=[gy gy gy gy gy gy];
    X = zeros(N,p);
    for(i=1:N)
        sigmay = Y(i);
        eps = mvnrnd(zeros(p+2,1),eye(p+2),1);
        X(i,:) = eps(1:p)' + sigmay*CS*eps((p+1):(p+2))';
    end
    X = X + 0.5*fY*beta*Gamma';
    
    for(h=3:hK)
        [WX,W]=fusing(Y,X,2,h);
        angles(h,isim) = subspace(W,CS)*180/pi;
    end
end



>> mean(angles')

ans =

  Columns 1 through 7

         0         0   36.8861   35.1913   34.5154   34.4559   34.3710

  Columns 8 through 14

   34.3689   34.3467   34.4344   34.4668   34.5156   34.5854   34.6636

  Column 15

   34.7515

>> std(angles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    0.5345    0.4774    0.4723    0.4750    0.4766

  Columns 8 through 14

    0.4847    0.4869    0.4869    0.4948    0.4937    0.4950    0.4981

  Column 15

    0.4983