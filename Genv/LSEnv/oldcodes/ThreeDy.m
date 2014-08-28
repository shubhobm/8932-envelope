function dims = ThreeDy(Xc,Yc,Sc,Sd,n,alpha)
p = size(Xc,2);
r = size(Yc,2);
dims = ones(3,1);
Sx = Sc(1:p,1:p);
Sy = Sc((p+1):(p+r),(p+1):(p+r));
Sxy = Sc(1:p,(p+1):(p+r));

newSc = [Sy, Sxy';Sxy, Sx];
newSd = blkdiag(Sy,Sx);
dims = ThreeDx(Yc,Xc,newSc,newSd,n,alpha);