function w = update_neural_network(w,data, radius)
%This function will take in neural node cooridinates w and update their
%position to be the average of the data around them
% We can model the neural network as network of springs each such that
% their relative distances in t0 is the x0 of the spring
step = 0.05;
for l = 1:100 % perform group average movement
    for i = 1:numel(w(:,1))
        % Determine error due to data position
        d_dist = data(:,1)*0;
        for k = 1:3
            d_dist = (data(:,k) - w(i,k)).^2 + d_dist;
        end
        index = d_dist.^0.5 < radius*2;
        if sum(index) > 1
            for k= 1:3
                dw(i,k) = mean(data(index,k)) - w(i,k); % This is the correction value based on the data
            end
        else
            dw(i,:) = [0, 0, 0];
        end
        
    end
    for k = 1:3

        w(:,k) = w(:,k) + step * mean(dw(:,k));
%         w(:,k) = w(:,k) + step/100 * dw(:,k);
    end

end
% for i = 1:10 % perform individual step
%         IDx = knnsearch(w,data);
%         for i = 1:numel(w(:,1))
%             index = IDx == i; % Grab all points whose NN is node i
% %             Determine error due to data position
%             d_dist = data(index,1)*0;
%             for k = 1:3
%                 d_dist = (data(index,k) - w(i,k)).^2 + d_dist;
%             end
%             d_index = d_dist.^0.5 > 2*radius;
%             index(d_index) = [];
%             if sum(index) > 20 
%                        
%             for k= 1:3
%                 if k == 3
%                     step = step/10;
%                 else
%                     step = 0.0;
%                 end
%                 dx(k) = mean(data(index,k)) - w(i,k); % This is the correction value based on the data
%                 w(i,k) = w(i,k) + step*dx(k);
%             end
%             end
%     
%         end
% end
    
end