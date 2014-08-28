%% Homework 2
load skulls.txt;
Y = skulls(:,1:4);
X = [(skulls(:,5)==-3300) ...
    (skulls(:,5)==-1850) ...
    (skulls(:,5)==-200) ...
    (skulls(:,5)==150)];

% Obtaining the best envelope dimension
aicmodel = modelselectaic(X,Y,'env')
bicmodel = modelselectbic(X,Y,'env')
lrtmodel = modelselectlrt(X,Y,0.05,'env')
% aic selects u=3, other two select u=1 so we go with u=1
u = 1;

% fitting envelope model
n = length(X)
env1 = env(X,Y,u); ols1 = fit_OLS(X,Y);
[env1.beta ols1.betaOLS]
%% First 4 cols give beta for envelope model while last 4 give beta for OLS
%     0.5633    2.5697    3.9983    5.0153    1.0000    3.1000    4.1333    4.8000
%    -0.2517   -1.1483   -1.7867   -2.2412   -0.9000    0.2000   -1.3000   -3.2667
%    -0.6584   -3.0031   -4.6727   -5.8612   -0.1000   -3.1333   -4.6333   -5.6667
%     0.1031    0.4705    0.7321    0.9183   -0.3000    0.0333    1.4333    0.8333
[env1.Sigma ols1.SigmaOLS]
%% First 4 cols give Sigma for envelope model while last for for OLS
%    22.9510    1.4750   -0.8249    3.3251   20.4071    0.0356    0.0764    1.9420
%     1.4750   21.7500    2.4735    2.5711    0.0356   22.7018    5.0267    2.7502
%    -0.8249    2.4735   22.1748   -0.2436    0.0764    5.0267   23.3731    1.0956
%     3.3251    2.5711   -0.2436   10.4252    1.9420    2.7502    1.0956    9.8142
env1.asySE
%     9.6409   10.7346   12.1860   13.4535
%     4.4661    7.1991    9.9708   12.0794
%    11.2326   11.8866   12.8008   13.6349
%     1.8951    3.7084    5.3778    6.6135
env1.covMatrix

env1.ratio
%     1.4817    1.3308    1.1723    1.0618
%     3.3737    2.0929    1.5111    1.2473
%     1.3611    1.2862    1.1943    1.1213
%     5.2275    2.6714    1.8422    1.4979

% doing hypothesis tests
stat1 = 30*(env1.beta(1,:))*inv(env1.covMatrix(1:4,1:4))*(env1.beta(1,:))'
pval1 = 1 - chi2cdf(stat1,145)
% stat1 =
% 
%   1.7247e+003
% 
% 
% pval1 =
% 
%      0

muNew = env1.beta(1,:) - env1.beta(2,:)
covNew = env1.covMatrix(1:4,1:4) + env1.covMatrix(5:8,5:8)
            - 2*env1.covMatrix(1:4,5:8)
stat2 = 60*muNew*inv(covNew)*muNew'
pval2 = 1 - chi2cdf(stat2,145)
% stat2 =
% 
%   332.8566
% 
% 
% pval2 =
% 
%   1.1102e-016