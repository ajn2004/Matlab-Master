function [locations, num_mol, a3] = neural_calc(i1, Theta1, Theta2)
% A function to calculate the probable location of molecules in an image i1
% given the neural net parameters theta1 and theta2
% AJN 10/23/15

[m, n] = size(i1); % calculate size of image

[Xgrid, Ygrid] = meshgrid(-3:3,-3:3);
locations = [];
a3 = zeros(m,n);
i2 = a3;
num_mol = 0;
locations = [];
for i = 4:m-3
    for j = 4:n-3
        X = double(i1(i-3:i+3,j-3:j+3));
        x = [ 1, X(:).'];
        a2 = sigmoid(x*Theta1.');
        y = [ones(numel(a2(:,1)),1), a2];
        a3(i,j) = sigmoid(y*Theta2.');
        xcm = sum(sum(X.*Xgrid))./sum(sum(X));
        ycm = sum(sum(X.*Ygrid))./sum(sum(X));
        %         if a3(i,j) > 0.1 && round(ycm) == 0 && round(xcm) == 0
        
    end
end

% a3 = medfilt2(a3, [1,1]);

for i = 4:m-3
    for j = 4:n-3
        if (a3(i,j) > 2.5*std(a3(:))+mean(a3(:))) && (a3(i,j) > a3(i-1,j-1)) && (a3(i,j) > a3(i-1,j)) && (a3(i,j) > a3(i-1,j+1)) && (a3(i,j) > a3(i,j-1)) && (a3(i,j) > a3(i,j+1)) && (a3(i,j) > a3(i+1,j-1)) && (a3(i,j) > a3(i+1,j)) && (a3(i,j) > a3(i+1,j+1))
            num_mol = num_mol +1;
            locations(num_mol,:) = [j,i];
        end
    end
end
