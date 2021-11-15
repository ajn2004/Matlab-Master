function upsampled_surface = upsample_smooth_surface(smooth_surface,kk,shrink, iterations)

[k, index] = knnsearch(smooth_surface,smooth_surface,'k',kk);
for j = 1:iterations
    upsampled_surface = [];
    
    for i = 1:numel(k(:,1))
        for l = 2:numel(k(1,:))
            upsampled_surface = [upsampled_surface; mean(smooth_surface(k(i,[1,l]),:))];
        end
        upsampled_surface = [upsampled_surface; mean(smooth_surface(k(i,2:end),:))];
    end
    
    k = boundary(upsampled_surface,shrink);
    for i = 1:numel(k(:,1))
        upsampled_surface = [upsampled_surface; mean(upsampled_surface(k(i,:),:))];
    end
    smooth_surface = [smooth_surface; upsampled_surface];    
    [k, index] = knnsearch(smooth_surface,smooth_surface,'k',kk);
%     for i = 1:numel(k(:,1))
%         upsampled_surface = [upsampled_surface; mean(smooth_surface(k(i,[1,6:10]),:))];
%     end
end
upsampled_surface = smooth_surface;

% Boundary Averaging Behavior
% function upsampled_surface = upsample_smooth_surface(smooth_surface,shrink, iterations)
% 'old method' adds a point at the center of every boundary triangle

% k = boundary(smooth_surface,shrink);
% for j = 1:iterations
%     upsampled_surface = [];
%     
%     for i = 1:numel(k(:,1))
%         upsampled_surface = [upsampled_surface; mean(smooth_surface(k(i,:),:))];
%     end
%     smooth_surface = [smooth_surface; upsampled_surface];
%     k = boundary(smooth_surface,shrink);
% end
% upsampled_surface = smooth_surface;
