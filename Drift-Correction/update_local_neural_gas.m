function w = update_local_neural_gas(w,data, radius)
%This function will take in neural node cooridinates w and update their
%position to be the average of the data around them
% We can model the neural network as network of springs each such that
% their relative distances in t0 is the x0 of the spring
step = 0.05;

for i = 1:10 % perform individual step
        IDx = knnsearch(w,data);
        for i = 1:numel(w(:,1))
            index = IDx == i; % Grab all points whose NN is node i
%             Determine error due to data position
            d_dist = data(index,1)*0;
            for k = 1:3
                d_dist = (data(index,k) - w(i,k)).^2 + d_dist;
            end
            d_index = d_dist.^0.5 > 2*radius;
            index(d_index) = [];
            if sum(index) > 20 
                       
            for k= 1:3
                
                dx(k) = mean(data(index,k)) - w(i,k); % This is the correction value based on the data
                w(i,k) = w(i,k) + step/10*dx(k);
            end
            end
    
        end
end
    
end