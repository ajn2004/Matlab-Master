% matlab nn player
close all; clearvars; clc;

net = feedforwardnet([20]);
net.Inputs{1}.size = 2;
view(net)