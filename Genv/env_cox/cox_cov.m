function [M U] = cox_cov(Y,X,delta,b)
N = size(X,1);
p = size(X,2);
wts = zeros(N,1);
% for i=1:N
%     vi = 0;
%     for j=1:N
%         if Y(j)>=Y(i)
%             vi = vi + exp(X(j,:)*b);
%         end
%     end
%     wts(i) = exp(X(i,:)*b)*delta(i)./vi;
% end

for i=1:N
    vi = (Y>=Y(i))'*exp(X*b);
    wts(i) = exp(X(i,:)*b)*delta(i)./vi;
end


wts = wts./mean(wts);
Exw = wts'*X/N;
Sxw = (X-Exw(ones(N,1),:))'*diag(wts)*(X-Exw(ones(N,1),:))./N;
f = X*b;
for i=1:N
    if wts(i)~=0
        f = f + (delta(i) - wts(i))/wts(i);
    end
end
Ys = f;
Eyw = wts'*Ys/N;
Sxyw = (X-Exw(ones(N,1),:))'*diag(wts)*(Ys-Eyw(ones(N,1),:))./N;
Syw = (Ys-Eyw(ones(N,1),:))'*diag(wts)*(Ys-Eyw(ones(N,1),:))./N;
M = Sxw;
U = Sxyw*Sxyw'./Syw;