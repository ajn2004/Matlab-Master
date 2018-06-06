function [ locations, num_mol, a3] = neural_gpu(a3)
% Support software to run the gpu neural net called image_neural_2
% AJN 11/4/15

[m, n, o] = size(a3); % calculate size of image
% [m, n, o] = size(a3); % calculate size of image



% a3 = zeros(m,n,o);

num_mol = 0;
locations = [];

% [a3] = image_neural_2(double(a3), Theta1, Theta2, o);

% a3 = medfilt2(a3, [1,1]);

% tic
for i = 4:m-3
    for j = 4:n-3
        if (a3(i,j) > 0.3) && (a3(i,j) > a3(i-1,j-1)) && (a3(i,j) > a3(i-1,j)) && (a3(i,j) > a3(i-1,j+1)) && (a3(i,j) > a3(i,j-1)) && (a3(i,j) > a3(i,j+1)) && (a3(i,j) > a3(i+1,j-1)) && (a3(i,j) > a3(i+1,j)) && (a3(i,j) > a3(i+1,j+1)) && a3(i,j) > a3(i,j-2) && a3(i,j) > a3(i,j-2)
            if a3(i,j) > a3(i-2,j-2) && a3(i,j) > a3(i-2,j-1) && a3(i,j) > a3(i-2,j) && a3(i,j) > a3(i-2,j+1) && a3(i,j) > a3(i-2,j+2) && a3(i,j) > a3(i-1,j-2) && a3(i,j) > a3(i-1,j+2) && a3(i,j) > a3(i-1,j+2)
                if a3(i,j) > a3(i+2,j-2) && a3(i,j) > a3(i+2,j-1) && a3(i,j) > a3(i+2,j) && a3(i,j) > a3(i+2,j+1) && a3(i,j) > a3(i+2,j+2) && a3(i,j) > a3(i+1,j-2) && a3(i,j) > a3(i+1,j+2) && a3(i,j) > a3(i+1,j+2)
                    num_mol = num_mol +1;
                    locations(num_mol,:) = [j,i];
                end
            end
        end
    end
end
% toc

