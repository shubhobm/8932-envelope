function [alpha Gamma Eta] = Logistic_ini1D(Y,X)
% [N p] = size(X);
% Sx = cov(X);
% ab = glmfit(X,Y,'binomial','link','logit','constant','on');
% a = ab(1);
% b = ab(2:end);
% [M U] = Logistic_cov(Y,X,a,b);
% [v1,d1] = eig(M);
% [v2,d2] = eig(M-U);
% v = [v1 v2];
% 
% W0 = v(:,1);
% theta = a + X*W0*inv(W0'*M*W0)*W0'*U;
% Fw0 = - Y'*theta + sum(log(1+exp(theta)))+ N/2*log(W0'*Sx*W0) + N/2*log(W0'*inv(Sx)*W0);
% 
% for i=2:(2*p)
%     W = v(:,i);
%     theta = a + X*W*inv(W'*M*W)*W'*U;
%     Fw = - Y'*theta + sum(log(1+exp(theta)))+ N/2*log(W'*Sx*W) + N/2*log(W'*inv(Sx)*W);
%     if Fw<Fw0
%         W0 = W;
%         Fw0 = Fw;
%     end
% end
% 
% ab = glmfit(X*W0,Y,'binomial','link','logit','constant','on');
% alpha = ab(1);
% Eta = ab(2);
% Gamma = W0;

ab = glmfit(X,Y,'binomial','link','logit','constant','on');
alpha = ab(1);
Gamma = orth(ab(2:end));
Eta = Gamma'*ab(2:end);
