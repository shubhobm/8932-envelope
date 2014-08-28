%% Homework 1
%% Q1
Sig = [1 0 1/2; 0 1 1/2; 1/2 1/2 1]
[vecs vals] = eig(Sig)

% Output

% Sig =
% 
%     1.0000         0    0.5000
%          0    1.0000    0.5000
%     0.5000    0.5000    1.0000
% 
% 
% vecs =
% 
%    -0.5000    0.7071   -0.5000
%    -0.5000   -0.7071   -0.5000
%     0.7071   -0.0000   -0.7071
% 
% 
% vals =
% 
%     0.2929         0         0
%          0    1.0000         0
%          0         0    1.7071

%% Q3
load MPLS46.txt
Y = MPLS46(:,1:4); n = size(MPLS46,1);
env2 = envmean(Y, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 1.3 on page 13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean(Y)'...
sqrt(diag(cov(Y)))./sqrt(n)...
env2.mu...
env2.asySE./sqrt(n)...
bstrp_envmean(Y, u, 100)...
sqrt(diag(cov(Y)))./env2.asySE]

% Basis matrices for envelope and its orthogonal compliment
env2.Gamma
env2.Gamma0

% Output:
% ans =
% 
%     0.2807    0.6760
%     0.6794   -0.4011
%     0.5747   -0.2179
%     0.3596    0.5785
% 
% 
% ans =
% 
%    -0.1127   -0.6720
%    -0.6142   -0.0167
%     0.7811   -0.1101
%          0    0.7322

env1 = envmean(Y,1);
env1.Gamma
env1.Gamma0

% Output:
% ans =
% 
%     0.4866
%     0.5105
%     0.4709
%     0.5300
% 
% 
% ans =
% 
%    -0.7238   -0.3831   -0.3041
%     0.6900   -0.4019   -0.3190
%          0    0.8317   -0.2943
%          0         0    0.8480