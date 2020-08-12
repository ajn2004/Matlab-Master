function plot_color_data(cdata)
% This is a quick function to plot orange and red data without an error if
% one or the other doesn't exist.
try
   plot3(cdata.red.xf,cdata.red.yf,cdata.red.zf,'.r')
   xlim([0 180]);
   ylim([0 180]);
   zlim([-1 1]/0.133)
   hold on
catch lsterr
    
end

try
    plot(cdata.orange.xf,cdata.orange.yf,cdata.orange.zf,'.b')
    hold off
catch lsterr
    
end
