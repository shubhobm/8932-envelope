function [dx,dy] = aict(Yc,Xc,mindim,Sc,Sd,n)
dx = mindim;
dy = mindim;
p = size(Xc,2);
r = size(Yc,2);
[Ghat,Hhat] = AA_env(Yc,Xc,dy,dx,Sc,Sd);
aic0 = p*(p+3)+r*(r+3)+2*dx*dy+2*logL(Sc,Ghat,Hhat,n);
        
        
for l = mindim:(p-1)
    for m = mindim:(r-1)
        [Ghat,Hhat] = AA_env(Yc,Xc,m,l,Sc,Sd);
        aicc = p*(p+3)+r*(r+3)+2*l*m+2*logL(Sc,Ghat,Hhat,n);
        if aicc < aic0
            aic0 = aicc;
            dx = l;
            dy = m;
        end
    end
end
