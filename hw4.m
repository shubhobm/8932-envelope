load AISxenv.txt
Y = AISxenv(:,1);
X = AISxenv(:,2:3);
[~,idx] = max(Y); Y(idx) = []; X(idx,:) = [];
maxdata = AISxenv(idx,:)

% select envelope dimension
u = modelselectaic(X, Y, 'xenv')
u = modelselectbic(X, Y, 'xenv')
alpha = 0.01;
u = modelselectlrt(X, Y, alpha,'xenv')

% Select u=1 since 2 is the full predictor dimension
u = 1; n = length(Y);

% Doing OLS and envelope fits, comparing models
envX = xenv(X, Y, u); olsX = fit_OLS(X,Y);
[olsX.betaOLS' envX.beta envX.asySE sqrt(olsX.n)*envX.beta./envX.asySE envX.ratio envX.Gamma]

% Checking excluded point

MSE = sum((Y - X*envX.beta).^2)/n
maxpred = (maxdata(1)-maxdata(2:3)*envX.beta)^2