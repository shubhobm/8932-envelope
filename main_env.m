% Applying envelope model with 0/1 as predictor and
% Other data as response

%% Preparing the data
load data508.txt;
names = textread('vars.txt','%s\t');
[len wid] = size(data508);
types = data508(1,2:wid-1);
X = data508(2:len,2:wid-1);
X1 = transform(X);

Xts = transform(X(:,find(types==1))); Xtc = transform(X(:,find(types==2)));
X3d = transform(X(:,find(types==3))); Xqc = transform(X(:,find(types==4)));

% Xts = X(:,find(types==1)); Xtc = X(:,find(types==2));
% X3d = X(:,find(types==3)); Xqc = X(:,find(types==4));

Y = data508(2:len,wid);
names = names(2:wid-1);

%% Applying envelope assuming 0/1 mutagenicity as predictor
mod_ts = env_pca(Xts,Y,90);
mod_tc = env_pca(Xtc,Y,90);
mod_tstc = env_pca([Xts Xtc],Y,90);

mod_all = env_pca(X1,Y,90);
% mod_allseq = envseq_pca(X1,Y,90);
% Analyze(Xts, Y, mod_ts, names(find(types==1)), types(find(types==1)),0); pause
% Analyze(Xtc, Y, mod_tc, names(find(types==2)), types(find(types==2)),0); pause
% Analyze([Xts Xtc], Y, mod_tstc, names(find(types==1|types==2)),
% types(find(types==1|types==2)),0); pause
% Analyze(X1, Y, mod_all, names, types,0)

%%
% pcX = X1*mod_allseq.loadings;
% xlin = min(pcX(:,1)):.01:max(pcX(:,1));
% g = mod_allseq.envMod.Gamma(1:2);
% ylin = xlin*g(1)/g(2);
% hold on
% gscatter(pcX(:,1),pcX(:,2),Y);
% plot(xlin,ylin);
% hold off
