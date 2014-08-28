setwd("C:/Study/UMN files/8932")
fulldata = read.table('data508.txt', header=F)
atom = read.table('data508atom1.txt', header=F)

require(spls); require(plsRglm)
d = dim(atom)
X = atom[2:d[1], 2:d[2]-1]
Y = atom[2:d[1], d[2]]
types = atom[1, 2:d[2]-1]

omitvars = c()
for(i in 1:(d[2]-1)){
  if(var(X[,i])==0){omitvars = c(omitvars,i)}
}
X = X[,-omitvars]

cv.spls(X[,500],Y,fold=10, eta = seq(0.1,0.9,0.1), K = c(1:10))
cv.sgpls(X,Y,fold=10, eta = seq(0.1,0.9,0.1), K = c(1:10))