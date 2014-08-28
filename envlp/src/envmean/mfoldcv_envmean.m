%% mfoldcv_meanenv
% Select the dimension of the envelope subspace using m-fold cross validation.

%% Syntax
%         u = mfoldcv_meanenv(Y, m)
%         u = mfoldcv_meanenv(Y, m, Opts)

%% Input
%
% *Y*: Data matrix. An n by p matrix, p is the dimension of the variable
% and n is number of observations. 
%
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.

%% Output
% 
%  *u*: The dimension of the envelope subspace selected by m-fold cross
%  validation.  An integer between 0 and p.

%% Description
% This function implements m-fold cross validation to select the dimension
% of the envelope space, based on prediction performance.  For each u, the
% data is partitioned into m parts, each part is in turn used for testing 
% for the prediction performance while the rest m-1 parts are used for 
% training.  The dimension is selected as the one that minimizes the average 
% prediction errors. As Y is multivariate, the identity inner product is 
% used for computing the prediction errors.

%% Example
% 
%         load Adopted
%         Y = Adopted(:, 1 : 6);
%         u = mfoldcv_envmean(Y, 5)


function u = mfoldcv_envmean(Y, m, Opts)

if nargin < 2
    error('Inputs: Y and m should be specified!');
elseif nargin == 2
    Opts = [];
end

Y = double(Y);

[n, r]=size(Y);

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

tempInd = min(floor((m - 1) * n / m) - 1, r);
PreErr = zeros(m, tempInd + 1);

for j = 0 : tempInd
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(j + 1) '\n']);
    end
    
    for i = 1 : m

        index = logical([ones(n, 1)]);
        index((floor((i - 1) * n / m) + 1) : ceil(i * n / m)) = 0;
        tempY = Y(index, :);
        ModelTemp = envmean(tempY, j);
        
        testY = Y(logical(1 - index), :);
        testN = size(testY, 1);
        resi = mean(testY) - ModelTemp.mu';
        PreErr(i, j + 1) = sqrt(resi * resi');
        
    end
end


[minErr ind] = min(mean(PreErr));
u = ind - 1;
