function [M U] = Logistic_cov(Y,X,a,b)
N = size(X,1);
p = size(X,2);
wts = 1./(2+exp(-a-X*b)+exp(a+X*b));
wts = wts./mean(wts);
Exw = wts'*X/N;
Sxw = (X-Exw(ones(N,1),:))'*diag(wts)*(X-Exw(ones(N,1),:))./N;
theta = a + X*b;
Ys = theta;
for i=1:N
    if wts(i)~=0
        Ys(i) = Ys(i) + (Y(i)-1./(1+exp(-theta(i))))./wts(i);
    end
end

Eyw = wts'*Ys/N;
Sxyw = (X-Exw(ones(N,1),:))'*diag(wts)*(Ys-Eyw(ones(N,1),:))./N;
Syw = (Ys-Eyw(ones(N,1),:))'*diag(wts)*(Ys-Eyw(ones(N,1),:))./N;
M = Sxw;
U = Sxyw*Sxyw'./Syw;