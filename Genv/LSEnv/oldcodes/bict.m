function [dx,dy] = bict(Yc,Xc,mindim,Sc,Sd,n)
dx = mindim;
dy = mindim;
p = size(Xc,2);
r = size(Yc,2);
[Ghat,Hhat] = AA_env(Yc,Xc,dy,dx,Sc,Sd);
bic0 = (p*(p+3)+r*(r+3)+2*dx*dy)/2*log(n)+2*logL(Sc,Ghat,Hhat,n);
bicttt = zeros(p,r);        
bicttt(1,1)=bic0;        
for l = mindim:(p-1)
    for m = mindim:(r-1)
        [Ghat,Hhat] = AA_env(Yc,Xc,m,l,Sc,Sd);
        bbic = (p*(p+3)+r*(r+3)+2*l*m)/2*log(n)+2*logL(Sc,Ghat,Hhat,n);
        bicttt(l+1,m+1) = bbic;
        if bbic < bic0
            bic0 = bbic;
            dx = l;
            dy = m;
        end
    end
end
