function get_coords_in_workspace(cdata)
% Get mean normalized and transformed coordinates into the matlab workspace
assignin('base','xfr',cdata.red.xf - mean(cdata.red.xf));
assignin('base','yfr',cdata.red.yf - mean(cdata.red.yf));
assignin('base','zfr',cdata.red.zf - mean(cdata.red.zf));
try
if min(cdata.orange.xf) > max(cdata.red.xf) % if the most left orange is more right than the furthest right red, apply transform
[xf, yf] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
else
    % otherwise, assume coords are already transformed and just assign
    xf = cdata.orange.xf;
    yf = cdata.orange.yf;
end
assignin('base','xfo',xf - mean(cdata.red.xf));
assignin('base','yfo',yf - mean(cdata.red.yf));
assignin('base','zfo',cdata.orange.zf - mean(cdata.red.zf));
catch
end
end