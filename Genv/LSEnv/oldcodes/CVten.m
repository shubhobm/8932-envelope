function [dx,dy] = CVten(Yc,Xc,dys,dxs)
p = size(Xc,2);
r = size(Yc,2);
N = size(Yc,1);
Ntest = floor(N/10);
Ntrain = N-Ntest;
Nfold = 10;

datavec = [Xc,Yc];

devs = zeros(length(dxs),length(dys));

        
for isim = 1:10
    data_train = datavec(:,:);
    idx = [1:Ntest] + (isim-1)*Ntest;
    data_train(idx,:)=[];
    data_test = datavec(idx,:);

    mu = mean(data_train);
    data_train = data_train - mu(ones(Ntrain,1),:);
    data_test = data_test - mu(ones(Ntest,1),:);

    Xc = data_train(:,1:p);
    Yc = data_train(:,(p+1):(p+r));

    XcN = data_test(:,1:p);
    YcN = data_test(:,(p+1):(p+r));

    Sc = data_train'*data_train./Ntrain;
    Sx = Sc(1:p,1:p);
    Sy = Sc((p+1):(p+r),(p+1):(p+r));
    Sxy = Sc(1:p,(p+1):(p+r));
    Sd = blkdiag(Sx,Sy);
    Soff = Sc - Sd;
    
    for dx = dxs
        for dy = dys
            [Ghat,Hhat] = AA_env(Yc,Xc,dy,dx,Sc,Sd);
            Ghat = orth(Ghat);
            Hhat = orth(Hhat);
            betahat = Ghat*inv(Ghat'*Sx*Ghat)*Ghat'*Sxy*Hhat*Hhat';
            resi = YcN - XcN*betahat;
            devs(dx,dy) = devs(dx,dy)+sqrt(trace(resi'*resi*resi'*resi));
            end
    end
end

i =1;
j =1;
dev0 = devs(1,1);

for i=dxs
    for j=dys
        if dev0 > devs(i,j)
            dx = dxs(i);
            dy = dys(j);
            dev0 = devs(i,j);
        end
    end
end

