%% bstrp_xenv
% Compute bootstrap standard error of the envelope model for the reduction on X. 

%% Syntax
%         bootse = bstrp_xenv(X, Y, u, B)
%         bootse = bstrp_xenv(X, Y, u, B, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations. The predictors must be continuous variables.
% 
% *Y*: Responses. An n by r matrix, r is the number of
% responses. The response can be univariate or multivariate and must be
% continuous variable.
% 
% *u*: Dimension of the envelope subspace.  A positive integer between 0
% and p.
% 
% *B*: Number of bootstrap samples.  A positive integer.
% 
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out the number of bootstrap samples, 
% logical 0 or 1. Default value: 0.
%
%% Output
%
% *bootse*: The standard error for elements in $$\beta$ computed by
% bootstrap.  An p by r matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the envelope model by bootstrapping the residuals. The
% envelope model here is for the reduction on X.

%% Example
%
%         load wheatprotein.txt
%         X = wheatprotein(:, 1 : 6);
%         Y = wheatprotein(:, 7);
%         alpha = 0.01;
%         u = lrt_xenv(X, Y, alpha)
%         B = 100;
%         bootse = bstrp_xenv(X, Y, u, B)

function bootse = bstrp_xenv(X, Y, u, B, Opts)

if (nargin < 4)
    error('Inputs: X, Y, u and B should be specified!');
elseif (nargin == 4)
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

X = double(X);
Y = double(Y);

[n r] = size(Y);
p = size(X, 2);

ModelOutput = xenv(X, Y, u, Opts);

Yfit = ones(n, 1) * ModelOutput.mu' + X * ModelOutput.beta;
resi = Y - Yfit;

bootBeta = zeros(B, p * r);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of bootstrap sample ' int2str(i) '\n']);
    end
    
    bootresi = resi(randsample(1 : n, n, true), :);
    Yboot = Yfit + bootresi;
    temp = xenv(X, Yboot, u, Opts);
    bootBeta(i, :) = reshape(temp.beta, 1, p * r);
    
end

bootse = reshape(sqrt(diag(cov(bootBeta, 1))), p, r);

fprintf('\nIf convergence is not reached for a bootstrap sample, \nit is still used in computing bootse.\n')