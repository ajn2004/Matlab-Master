function [x_out, y_out] = make_nn_channel_transform(x_in,y_in)

comp_name = get_computer_name();
load([comp_name,'\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_neural_net.mat'],'net1');
try
load([comp_name,'\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\zcalib.mat'],'q');
catch
    q = 0.122;
end
output = net1([x_in/q,y_in/q]');
x_out = output(1,:).'*q;
y_out = output(2,:).'*q;