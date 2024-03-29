%% predict_env
% Perform estimation or prediction under the envelope model.

%% Syntax
%         PredictOutput = predict_env(ModelOutput, Xnew, infType)
%
%% Input
%
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherited from env.
% 
% *Xnew*: The value of X with which to estimate or predict Y.  A p by 1
% vector.
% 
% *infType*: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
% 
%% Output
%
% *PredictOutput*: A list containing the results of the inference.
%
% * PredictOutput.value: The fitted value or the prediction value evaluated at
% Xnew. An r by 1 vector.
% * PredictOutput.covMatrix: The covariance matrix of PredictOutput.value. An r by r
% matrix.
% * PredictOutput.SE: The standard error of elements in PredictOutput.value.  An r
% by 1 vector. 

%% Description
% This function evaluates the envelope model at new value Xnew.  It can
% perform estimation: find the fitted value when X = Xnew, or prediction:
% predict Y when X = Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1 : 6);
%         alpha = 0.01;
%         u = lrt_env(X, Y, alpha);
%         ModelOutput = env(X, Y, u);
%         Xnew = X(2, :)';
%         PredictOutput = predict_env(ModelOutput, Xnew, 'estimation')
%         [PredictOutput.value, Y(1, :)']  % Compare the fitted value with the data
%         PredictOutput.SE
%         PredictOutput = predict_env(ModelOutput, Xnew, 'prediction')
%         PredictOutput.value
%         PredictOutput.SE

function PredictOutput = predict_env(ModelOutput, Xnew, infType)

if nargin < 3
    error('Inputs: ModelOutput, Xnew and infType should be specified!');
end

if ~strcmp(infType, 'estimation') && ~strcmp(infType, 'prediction')
    error('Inference type can only be estimation or prediction.');
end

[r p] = size(ModelOutput.beta);
[s1 s2] = size(Xnew);

if s1 ~= p || s2 ~= 1
    error('Xnew must be a p by 1 vector');
end

n = ModelOutput.n;
u = size(ModelOutput.Gamma, 2);

if u == 0
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.alpha;
        PredictOutput.covMatrix = ModelOutput.Sigma / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.alpha;
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.Sigma;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
else
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.alpha + ModelOutput.beta * Xnew;
        PredictOutput.covMatrix = ModelOutput.Sigma / n ...
            + kron(Xnew', eye(r)) * ModelOutput.covMatrix * kron(Xnew, eye(r)) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.alpha + ModelOutput.beta * Xnew;
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.Sigma ...
            + kron(Xnew', eye(r)) * ModelOutput.covMatrix * kron(Xnew, eye(r)) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
end