function [G,H] = alter_env(Y,X,dimy,dimx,Sc,Sd,Go,Ho)
% G is the basis for the sigma_X envelope
% H is the basis for the sigma_Y envelope
% for fixed G, the objective function is the same as the objective function
% of envelope model of Y on G'*X
% for fixed H, the objective function is the same as the objective function
% of envelope model of X on H'*Y


G=Go;
H=Ho;

    [HX,H,L,dhat] = mlm_fit(Y,X*G,dimy);
    [GY,G,L,dhat] = mlm_fit(X,Y*H,dimx);