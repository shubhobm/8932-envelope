%% dF4henv
% The first derivative of the objective function for computing the envelope
% subspace in the heteroscedastic envelope model.

%% Syntax
%         df = dF4henv(R, DataParameter)
% 
%% Input
%
% *R*: An r by u semi orthogonal matrix, 0 < u <= r.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *df*: An r by u matrix containing the value of the derivative function evaluated at R.

%% Description
%
% The objective function is derived in Section 2.2 in Su and Cook (2012) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF4henv(R, DataParameter)

p = DataParameter.p;
n = DataParameter.n;
ng = DataParameter.ng;
sigRes = DataParameter.sigRes;
sigY = DataParameter.sigY;

df = zeros(size(R));

for i = 1 : p
    a = 2 * ng(i) / n * sigRes(:, :, i) * R * inv(R' * sigRes(:, :, i) * R);
    df = df + a;
end

b = 2 * inv(sigY) * R * inv(R' * inv(sigY) * R);

df = n * (df + b);