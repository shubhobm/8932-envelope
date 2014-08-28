function df = dF4cox_num(Y,Fh)
[N,P] = size(Y);
ep = 1e-6;
for k=1:P
    for j=1:N
        Yp = Y; Yp(j,k) = Yp(j,k)+ep;
        Ym = Y; Ym(j,k) = Ym(j,k)-ep;
        df(j,k) = (Fh(Yp)-Fh(Ym))/ep/2;
    end
end