function upsampled_surface = upsample_smooth_surface(smooth_surface,kk,shrink,m, iterations)

[k, index] = knnsearch(smooth_surface,smooth_surface,'k',kk);
for j = 1:iterations
    upsampled_surface = [];
    
%     for i = 1:numel(k(:,1))
%         for l = 2:numel(k(1,:))
%             upsampled_surface = [upsampled_surface; mean(smooth_surface(k(i,[1,l]),:))];
%         end
%         upsampled_surface = [upsampled_surface; mean(smooth_surface(k(i,2:end),:))];
%     end
%     upsampled_surface = unique(upsampled_surface,'rows');
%     smooth_surface = [smooth_surface; upsampled_surface]; 
%     for i = 1:3
%     smooth_surface(isnan(smooth_surface(:,i)),:) = [];
%     end
    k = boundary(smooth_surface,shrink);
    for i = 1:numel(k(:,1))
        smooth_surface = [smooth_surface; mean(smooth_surface(k(i,:),:))];
    end   
    [k, index] = knnsearch(smooth_surface,smooth_surface,'k',kk);

    % knn 'push'
%     ks = 3;
    [Idx, d] = knnsearch(smooth_surface,smooth_surface,"K",kk);
    push_surface = smooth_surface;

%     for i = 1:numel(Idx(:,1))
%         kick = [0,0,0];
%         for l = 1:3
%             for j = 2:kk
%                 kick(l) =kick(l)- m*diff(smooth_surface(Idx(i,[1,j]),l))./abs(diff(smooth_surface(Idx(i,[1,j]),l)))^1;
%             end
%         end
%         push_surface(i,:) = push_surface(i,:)+ kick/kk;
%     end
%     for i = 1:3
%         push_surface(isnan(push_surface(:,i)),:) = [];
%     end
%     smooth_surface = push_surface; 
%     k = boundary(smooth_surface,shrink);
%     for i = 1:numel(k(:,1))
%         smooth_surface = [smooth_surface; mean(smooth_surface(k(i,:),:))];
%     end
% %     smooth_surface = upsampled_surface;
end

k = boundary(smooth_surface,shrink);
uk = unique(k);



upsampled_surface = smooth_surface(uk,:);

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
