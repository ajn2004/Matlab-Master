% Trajec Viewer
close all
s = scatter3(xf_all*q,yf_all*q,zf_all*q);

axis equal
hold on
for i = 1:numel(trajec)
    ind = trajec(i).t;
    rgb = get_color(numel(ind));
    plot3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind));
    si = scatter3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind),[],rgb);
    si.MarkerFaceColor = si.MarkerEdgeColor;
end
axis equal

