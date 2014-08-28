function [dxs,dys] = tritests(Yc,Xc,mindim,Sc,Sd,n,alpha)
dx0 = mindim;
dy0 = mindim;
p = size(Xc,2);
r = size(Yc,2);
[Ghat,Hhat] = AA_env(Yc,Xc,dy0,dx0,Sc,Sd);

bic0 = (p*(p+3)+r*(r+3)+2*dx0*dy0)/2*log(n)+2*logL(Sc,Ghat,Hhat,n);
bic1 = zeros(p,r);        
bic1(1,1)=bic0;    

aic0 = p*(p+3)+r*(r+3)+2*dx0*dy0+2*logL(Sc,Ghat,Hhat,n);
aic1 = zeros(p,r);        
aic1(1,1)=aic0;    

Lik0 = 2*logL(Sc,Ghat,Hhat,n);
Lik1 = zeros(p,r);
Lik1(1,1) = Lik0;

for l = mindim:min(5,p-2)
    for m = mindim:min(5,r-2)
        [Ghat,Hhat] = AA_env(Yc,Xc,m,l,Sc,Sd);
        bic1(l+1,m+1) = (p*(p+3)+r*(r+3)+2*l*m)/2*log(n)+2*logL(Sc,Ghat,Hhat,n);
        aic1(l+1,m+1) = p*(p+3)+r*(r+3)+2*l*m+2*logL(Sc,Ghat,Hhat,n);
        Lik1(l+1,m+1) = 2*logL(Sc,Ghat,Hhat,n);
        fprintf('l=%d m=%d\n',l,m);
    end
end


for l = mindim:min(5,p-2)
    for m = mindim:min(5,r-2)
        if bic1(l+1,m+1) < bic0
            bic0 = bic1(l+1,m+1);
            dx_bic = l;
            dy_bic = m;
        end
        if aic1(l+1,m+1) < aic0
            aic0 = aic1(l+1,m+1);
            dx_aic = l;
            dy_aic = m;
        end
    end
end

dx_lrt = mindim;
dy_lrt = mindim;
for isim = 1:3
    l = dx_lrt;
    m = dy_lrt;
    for l=(dx_lrt+1):min(5,(p-2))
        chist = - Lik1(l+1,m+1) + Lik1(dx_lrt+1,m+1);
        df = m*(l-dx_lrt);
        if chist > chi2inv(1-alpha,df)
            dx_lrt = l;
            break;
        end
    end
    
    l = dx_lrt;
    m = dy_lrt;
    for m=(dy_lrt+1):min(5,r-2)
        chist = - Lik1(l+1,m+1) + Lik1(l+1,dy_lrt+1);
        df = l*(m-dy_lrt);
        if chist > chi2inv(1-alpha,df)
            dy_lrt = m;
            break;
        end
    end   
end

dxs = [dx_aic,dx_bic,dx_lrt];
dys = [dy_aic,dy_bic,dy_lrt];
