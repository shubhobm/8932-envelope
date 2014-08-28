load data508.txt;
[len wid] = size(data508);
types = data508(1,2:wid-1);
X = data508(2:len,2:wid-1);
X1 = transform(X);
Xts = transform(X(:,find(types==1))); Xtc = X(:,find(types==2));
X3d = X(:,find(types==3)); Xqc = transform(X(:,find(types==4)));
Y = data508(2:len,wid); Y = (Y-mean(Y))/std(Y);

%% Fitting envelopes
% k = xenv_pca(X1,Y,90);
% sigs = find(abs(k.tRatios)>2);
% length(sigs)
% [(types(sigs))' k.tRatios(sigs)]

%% Envelope predictions
testsize = ceil((len-1)/10);
tests = (vec2mat(randsample(RandStream.create('mt19937ar','NumStreams',1), len-1, len-1, false),testsize))';

cnt = 0; predY = zeros(len-1,2); dims = zeros(9,1);
for i = 1:9
    trnmod = xenv_pca(removerows(Xts,tests(:,i)),removerows(Y,tests(:,i)),95);
    predY((cnt+1):(cnt+testsize)) = Xts(tests(:,i),:)*trnmod.beta;
    dims(i) = trnmod.u;
    cnt = cnt+testsize;
    fprintf('Done!%d\n',i);
%         if predval>0
%             predY(cnt) = 1;
%         else
%             predY(cnt) = 0;
%         end
end
predY(:,2) = Y;
dims