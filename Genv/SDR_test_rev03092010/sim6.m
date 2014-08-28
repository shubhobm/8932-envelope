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
        sigmay = gY(i,:)*[1 2 5]';
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

         0         0   73.4974   66.0922   51.0373   61.5553   60.0862

  Columns 8 through 14

   59.7128   60.2306   60.5020   59.7061   60.1384   61.4510   61.9224

  Column 15

   62.7898

>> mean(fangles')

ans =

  Columns 1 through 7

         0         0   73.4974   68.3107   59.0808   57.1870   54.5149

  Columns 8 through 14

   55.3238   53.3381   52.9566   52.8234   51.9425   50.7329   51.6717

  Column 15

   51.0228


>> std(angles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    1.2464    1.5259    1.5220    1.6614    1.7129

  Columns 8 through 14

    1.6029    1.5901    1.6586    1.5155    1.6059    1.4204    1.4130

  Column 15

    1.4843

>> std(fangles')/sqrt(nsim)

ans =

  Columns 1 through 7

         0         0    1.2464    1.4620    1.6591    1.7501    1.6112

  Columns 8 through 14

    1.6533    1.5750    1.5793    1.5814    1.4805    1.4988    1.4843

  Column 15

    1.4446



