function z = get_nn_z(iloc)
load('C:\Users\AJN Lab\Dropbox\Codes\Matlab Testing Folder\Simulations\Neural Net Sims\Out_Thetas.mat');
% for i = 1:o
%     isub = iloc(:,:,i);
%     x(i,:) = isub(:).';
% end
x = iloc.';
As = NN_eval(x,theta1,theta2);


z = As * ysc + ymi;
