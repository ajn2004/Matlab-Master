%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Network training routine to identify psfs
% By using a neural network we may be able to determine with some accuracy
% whether or not an image is a PSF
%
%
%
% AJN 10/7/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

[fname, fpath]= uigetfile('*.mat');
current_path = pwd;
addpath(pwd);
cd(fpath);
f_info = dir('*set.mat');
trainper = .6;


input_layer = 49;
hidden_layer = 100;

output = 1;
epsint = 0.12;


for i = 1:numel(f_info)
    load(f_info(i).name);
    if i == 1
        X = double(x); %prep data set for analysis
        Y = y;
    else
        X = vertcat(X,double(x));
        Y =vertcat(Y,y);
    end
end

load('C:\Users\FPALM8\Desktop\Matlab Testing Folder\PSF Learner\Simulator\Training_set.mat');
X = [X ; double(i3)];
Y = [Y; y];

        int_Theta1 = rand(hidden_layer, input_layer + 1)*2*epsint - epsint;  %initialize parameter set as 0s
        int_Theta2 = rand(output, hidden_layer +1)*2*epsint - epsint;
        lambda = 3;   %regulation factor, play around to find something useful based on trials not apart of the example set
        int_thetas = [int_Theta1(:); int_Theta2(:)];
        
        [J , grad] = neuralcostfunc(int_thetas, input_layer, hidden_layer, output, X, Y, lambda); %cost function to show first and final cost
        split = round(trainper * numel(X(:,1)));
        % conv = 0.001; %convergence factor, ask J to be within 100*conv % of itself between iterations before breaking analysis
        % [fin_thetas, i, Jout] = newtraphs(int_thetas, X, y.', lambda, conv); %apply netwon raphson method to determine proper values for theta parameters
        randind = randperm(numel(X(:,1)));
        Xtrain = X(randind(1:split),:);
        Ytrain = Y(randind(1:split));
        options = optimset('MaxIter', 50);
        [thetas] = fmincg(@(t)(neuralcostfunc(t, input_layer, hidden_layer, output, Xtrain, Ytrain, lambda)), int_thetas, options);
        
        Theta1 = reshape(thetas(1: hidden_layer*(input_layer +1)), hidden_layer, input_layer + 1);
        Theta2 = reshape(thetas(hidden_layer*(input_layer +1) +1: end), output, hidden_layer + 1);
        
        %% Determine accuracy
        X1 = [ ones(numel(X(:,1)),1), X];
        z2 = X1*Theta1.';
        a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
        z3 = a2*Theta2.';
        a3 = sigmoid(z3);
        Xtest = X(randind(split+1:end),:);
        Ytest = Y(randind(split+1:end));
        mean(a3 > 0.5 == Y) %

% save([current_path, 'neural_theta.mat'],'Theta1','Theta2');
save([current_path,'Neural_thetas_2.mat'],'Theta1','Theta2');