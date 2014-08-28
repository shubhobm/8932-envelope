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
fangles = zeros(hK,nsim);

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
        sigmay = fY(i,:)*[1 2 3 4 5]';
        eps = mvnrnd(zeros(p+2,1),eye(p+2),1);
        X(i,:) = eps(1:p)' + sigmay*CS*eps((p+1):(p+2))';
    end
    X = X + 0.5*fY*beta*Gamma';
    
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

         0         0   46.9682   46.1159   43.9323   45.5006   46.1129

  Columns 8 through 14

   46.9114   46.9702   47.3585   48.4249   49.4326   49.8974   52.2634

  Column 15

   53.4526

>> mean(fangles')

ans =

  Columns 1 through 7

         0         0   46.9682   45.1022   43.7012   43.6060   43.6690

  Columns 8 through 14

   43.7849   43.7502   43.8128   43.8976   43.9803   44.0411   44.1175

  Column 15

   44.2691

>> std(angles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    0.6009    0.5960    0.5692    0.5304    0.5375

  Columns 8 through 14

    0.6405    0.6675    0.6612    0.7043    0.6840    0.6802    0.8209

  Column 15

    0.8899

>> std(fangles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    0.6009    0.5743    0.5536    0.5205    0.5051

  Columns 8 through 14

    0.5141    0.5192    0.5267    0.5192    0.5152    0.5125    0.5149

  Column 15

    0.5161