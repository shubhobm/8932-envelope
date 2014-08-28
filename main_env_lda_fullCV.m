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
predY = zeros(len-1,1);
Xmod = [Xtc]; % change to any set of predictors

tic
textprogressbar('Doing leave-one-out CV: ');
for i = 1:len-1
    trnX = removerows(Xmod,i);
    trnY = removerows(Y,i);
    testX = Xmod(i,:);
    trnmod = env_pca(trnX,trnY,90,0);
    predY(i) = classify(testX*trnmod.loadings*trnmod.envMod.Gamma,...
    trnX*trnmod.loadings*trnmod.envMod.Gamma,trnY); 
    textprogressbar(i/(len-1)*100); % tracking loop progress
end
textprogressbar(' done');
toc

% predY(460:508) = xenv_pca_pred(removerows(Xmod,tests(1:49,10)),removerows(Y,tests(1:49,10)),95,Xmod(tests(1:49,10),:));
% fprintf('Done!10\n');
[ length(find(predY==1 & Y==1))/length(find(Y==1)) length(find(predY==0 & Y==0))/length(find(Y==0))]