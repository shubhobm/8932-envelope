function [pred] = xenv_pca_pred(X,Y,mod,Xnew)
% prediction by linear discriminant
X1 = X*mod.Gamma; Xn = Xnew*mod.Gamma;
fprintf('Done!1')
S0 = cov(X1(find(Y==0),:));
S1 = cov(X1(find(Y==1),:));
m0 = mean(X1(find(Y==0),:));
m1 = mean(X1(find(Y==1),:));
w = inv(S0+S1)*(m1-m0)';
c = w*(m0+m1)/2;

n = size(Xn,2);
pred = zeros(n,1);
for i = 1:n
    if w'*Xn(:,i) > c
        pred(i) = 1;
    else
        pred(i) = 0;
    end
end
