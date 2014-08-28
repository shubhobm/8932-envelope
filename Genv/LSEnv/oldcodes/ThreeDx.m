function dims = ThreeDx(Xc,Yc,Sc,Sd,n,alpha)
p = size(Xc,2);
r = size(Yc,2);
dims = ones(3,1);
bicval = zeros(p-2,1);
aicval = zeros(p-2,1);
likval = zeros(p-2,1);

for dn = 1:min(10,p-2)
    Ghat = X_env(Xc,Yc,Sc,dn);
    bicval(dn) = (p*(p+3)+r*(r+3)+2*dn*r)/2*log(n)+2*logL(Sc,Ghat,eye(r),n);
    aicval(dn) = p*(p+3)+r*(r+3)+2*dn*r+2*logL(Sc,Ghat,eye(r),n);
    likval(dn) = 2*logL(Sc,Ghat,eye(r),n);
end

bicval0 = bicval(1);
aicval0 = aicval(1);

for dn = 2:min(10,p-2)
    if bicval(dn) < bicval0
        bicval0 = bicval(dn);
        dims(2) = dn;
    end
    if aicval(dn) < aicval0
        aicval0 = aicval(dn);
        dims(1) = dn;
    end
end

for dn=2:min(10,p-2)
    chist = - likval(dn) + likval(dims(3));
    df = r*(dn-dims(3));
    if chist > chi2inv(1-alpha,df)
        dims(3) = dn;
    else
        break;
    end
end