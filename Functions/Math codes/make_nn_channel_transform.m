function [x_out, y_out] = make_nn_channel_transform(x_in,y_in)
comp_name = get_computer_name();
load([comp_name,'\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_neural_net.mat'],'net1');
output = net1([x_in,y_in]');
x_out = output(1,:).';
y_out = output(2,:).';