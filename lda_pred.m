function [pred] = lda_pred(X,Y,Xnew)
% prediction by linear discriminant
S0 = cov(X(find(Y==0),:));
S1 = cov(X(find(Y==1),:));
m0 = mean(X(find(Y==0),:));
m1 = mean(X(find(Y==1),:));
w = inv(S0+S1)*(m1-m0)';
c = w*(m0+m1)/2;

n = size(Xnew,1);
pred = zeros(n,3);
for i = 1:n
    if w*Xnew(i,:) > c
        pred(i,:) = [1 w*Xnew(i,:) c];
    else
        pred(i,:) = [0 w*Xnew(i,:) c];
    end
end