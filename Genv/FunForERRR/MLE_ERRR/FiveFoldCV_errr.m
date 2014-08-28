function [dev_ols dev_rrr dev_env dev_errr] = FiveFoldCV_errr(X,Y,d,u,maxitera,nsim)
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

Ntest = floor(N/5);
Ntrain = N-Ntest;
Nfold = 5;

datavec = [Xc,Yc];
dev_ols = 0;
dev_rrr = 0;
dev_env = 0;
dev_errr = 0;

for m = 1:nsim
    idx = randperm(N);
    datavec = datavec(idx,:);
    % this two lines permutes the original order of the dataset

    for isim = 1:Nfold
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

        %%%%%%%%%%%%%%%%%% OLS %%%%%%%%%%%%%%%%%%%%%%%
        betahat = Syx*inv(Sx);
        resi = YcN - XcN*betahat';
        dev_ols =  dev_ols + sqrt(trace(resi'*resi));

        %%%%%%%%%%%%%%%%%% Y-Env %%%%%%%%%%%%%%%%%%%%%%%
        [Ghat Ehat] = Yenv(Xc,Yc,u,maxitera);
        betahat = Ghat*Ehat;
        resi = YcN - XcN*betahat';
        dev_env =  dev_env + sqrt(trace(resi'*resi));

        %%%%%%%%%%%%%%%%%% RRR %%%%%%%%%%%%%%%%%%%%%%%
        [Ai Bi Aml Bml] = RRR(Xc,Yc,d);
        betahat = Aml*Bml';
        resi = YcN - XcN*betahat';
        dev_rrr = dev_rrr + sqrt(trace(resi'*resi));

        %%%%%%%%%%%%%%%%%% ERRR %%%%%%%%%%%%%%%%%%%%%%%
        [A Bi Ci Bml Cml] = YenvRRR(Xc,Yc,u,d,maxitera);
        betahat = A*Bml*Cml';
        resi = YcN - XcN*betahat';
        dev_errr = dev_errr  + sqrt(trace(resi'*resi));
    end
end

