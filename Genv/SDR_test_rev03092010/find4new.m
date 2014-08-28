function [] = find4new(varargin)
% This function reads additional inputs parameters for the EPFC
% model during call to function LDR.
% =============================================================
global NEWparameters;
NEWparameters.alpha = 0.05;
NEWparameters.fy = [];
NEWparameters.init = [];
for i=1:length(varargin),
    input = varargin{i};
    if ischar(input),
        switch lower(input),
            case 'alpha',
                NEWparameters.alpha = varargin{i+1};
            case 'fy',
                aux = varargin{i+1};
                for j=1:size(aux,2),
                    aux(:,j) = aux(:,j) - mean(aux(:,j));
                end
                NEWparameters.fy = aux;
            case 'initval',
               NEWparameters.init = varargin{i+1};
         end
    end
end



