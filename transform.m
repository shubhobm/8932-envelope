function [trmat] = transform(mat)
s = size(mat);
trmat = mat;
% fprintf('Start..')
for i = 1:s(1)
    for j = 1:s(2)
        if mat(i,j)<0
            C = -floor(mat(i,j));
        else
            C = 1;
        end
        trmat(i,j) = log(mat(i,j) + C);
    end
end

% for j = 1:s(2)
%     mj = mean(trmat(:,j));
%     sj = std(trmat(:,j));
%     for i = 1:s(1)
%         trmat(i,j) = (trmat(i,j)-mj)/sj;
%     end
% end
% fprintf('End.\n')