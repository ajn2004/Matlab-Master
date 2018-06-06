function [y, ymins, yscales] = func_feat_scale(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Scaling Function
%
% a short function to scale a data set so that the outputs are between 1
% and 0.
%
% ajn 12/22/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale columns
for i = 1:numel(y(1,:))
    ymins(1,i) = min(y(:,i));
    yscales(1,i) = (max(y(:,i)) - min(y(:,i)));
    y(:,i) = (y(:,i) - ymins(:,i))/yscales(:,i);
end
end