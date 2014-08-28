load data508.txt;
[len wid] = size(data508);
types = data508(1,2:wid-1);
X = data508(2:len,2:wid-1);
X1 = transform(X);

Xts = transform(X(:,find(types==1))); Xtc = transform(X(:,find(types==2)));
X3d = transform(X(:,find(types==3))); Xqc = transform(X(:,find(types==4)));

% Xts = X(:,find(types==1)); Xtc = X(:,find(types==2));
% X3d = X(:,find(types==3)); Xqc = X(:,find(types==4));

Y = data508(2:len,wid);

%% Envelope predictions
testsize = ceil((len-1)/10);
tests = (vec2mat(randsample(RandStream.create('mt19937ar','NumStreams',1), len-1, len-1, false),testsize))';
cnt = 0; predY = zeros(459,2);
Xmod = [Xtc];

for i = 1:9
    trnX = removerows(Xmod,tests(:,i));
    trnY = removerows(Y,tests(:,i));
    testX = Xmod(tests(:,i),:);
    trnmod = xenv_pca(trnX,trnY,90);
    predY((cnt+1):(cnt+testsize)) = classify(testX*trnmod.loadings*trnmod.envMod.Gamma,...
    trnX*trnmod.loadings*trnmod.envMod.Gamma,trnY); 
    cnt = cnt+testsize;
    fprintf('Done!%d\n',i);
end
% predY(460:508) = xenv_pca_pred(removerows(Xmod,tests(1:49,10)),removerows(Y,tests(1:49,10)),95,Xmod(tests(1:49,10),:));
% fprintf('Done!10\n');
predY(:,2) = Y(1:459);
length(find(predY(:,1)==1 & predY(:,2)==1))/length(find(predY(:,2)==1))
length(find(predY(:,1)==0 & predY(:,2)==0))/length(find(predY(:,2)==0))