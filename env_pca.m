function [outs] = env_pca(X,Y,pct,trace)

if nargin==3
    trace=1;
end

[a b c] = pcacov(X);
d = cumsum(c);
dim = min(find(d>pct));
loadA = a(:,1:dim);
trX = X*loadA;
% plot(trX(:,1),trX(:,2)); xlabel('PC1'); ylabel('PC2');
% outs.u = 10;
Opts.maxIter = 1000;
Opts.ftol = 1e-7;
Opts.gradtol = 1e-4;

if(trace==1)
    fprintf('Selecting envelope dimension.\n')
end
outs.u = modelselectlrt(Y,trX,0.05,'env',Opts);
if(trace==1)
    fprintf('Building envelope model..\n')
end
envmod = env(Y,trX,outs.u,Opts);

outs.envMod = envmod;
outs.pctvars = c(1:dim);
outs.loadings = loadA;
outs.beta = loadA*envmod.beta;
outs.asySE = sqrt((loadA.^2)*(envmod.asySE.^2));
outs.tRatios = outs.beta.*sqrt(length(X))./outs.asySE;
