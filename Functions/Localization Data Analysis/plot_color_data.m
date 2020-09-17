function plot_color_data(cdata)
% This is a quick function to plot orange and red data without an error if
% one or the other doesn't exist.
try
   % Create a graph to visualize red and Orange data
   ind_r = abs(cdata.red.zf) < 0.5;
   ind_o = abs(cdata.orange.zf) < 0.5;
   plot3(cdata.red.xf(ind_r),cdata.red.yf(ind_r),cdata.red.zf(ind_r),'.')
   hold on
   vector = xy_feature(cdata.orange.xf/cal.q,cdata.orange.yf/cal.q);
   xf_o = cal.q*cal.o2rx.'*vector.';
   yf_o = cal.q*cal.o2ry.'*vector.';
   plot3(xf_o(ind_o), yf_o(ind_o),cdata.orange.zf(ind_o),'.')
   hold off
   axis equal
   title('Non-corrected overlayed 3d localization data')
   legend('Glut4','vGlut1')
   correction = cdata.orange.zf - cdata.orange.zf_raw;
   plot3(xf_o(ind_o), yf_o(ind_o),cdata.orange.zf(ind_o) - correction(ind_o)*l,'.')
   % plot(xf_o(ind_o), yf_o(ind_o),'.')
   hold off
   axis equal
   % title('Fixed BEAS w/ tubulin ')
   legend('PA-JF','Glitter Bomb','Location',"northwest")
   % legend('Glitter Bomb','Location',"northwest")
   xlabel('Microns')
   ylabel('Microns')
catch lsterr
    
end

try
    plot(cdata.orange.xf,cdata.orange.yf,cdata.orange.zf,'.b')
    hold off
catch lsterr
    
end
