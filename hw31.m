%% HW3 Q1
load fiberpaper.dat
Y = fiberpaper(:, 1 : 4);
X = fiberpaper(:, [7 5 6]);
alpha = 0.01;
X1 = X(1, :)'; X4 = X(4, :)';

X_1 = X(setdiff(1:size(X,1),1),:); Y_1 = Y(setdiff(1:size(Y,1),1),:);
u1 = lrt_env(X_1, Y_1, alpha);
ModelOutput1 = env(X_1,Y_1, u1);
p1env = predict_env(ModelOutput1, X1, 'prediction');
p1penv = predict_env2(X_1, Y_1, X1, 'prediction');
[p1env.value p1env.SE p1penv.value p1penv.SE p1env.SE./p1penv.SE]

X_4 = X(setdiff(1:size(X,1),4),:); Y_4 = Y(setdiff(1:size(Y,1),4),:);
u4 = lrt_env(X_4, Y_4, alpha);
ModelOutput4 = env(X_4,Y_4, u4);
p4env = predict_env(ModelOutput4, X4, 'prediction');
p4penv = predict_env2(X_4, Y_4, X4, 'prediction');
[p4env.value p4env.SE p4penv.value p4penv.SE p4env.SE./p4penv.SE]

% OUTPUT
% 
% u =
% 
%      2
% 
% 
% ans =
% 
%    21.0006    2.6388   21.0352    2.6254    1.0051
%     7.0854    0.7162    7.0385    0.6966    1.0282
%     5.3011    1.2500    5.2725    1.2356    1.0116
%     0.8613    0.5833    0.8633    0.5711    1.0212
% 
% 
% u =
% 
%      2
% 
% 
% ans =
% 
%    21.8770    2.5521   21.8708    2.5608    0.9966
%     7.3094    0.6952    7.3202    0.6792    1.0236
%     5.7082    1.2082    5.7123    1.2017    1.0054
%     1.0508    0.5639    1.0571    0.5584    1.0099


