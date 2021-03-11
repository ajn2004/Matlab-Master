function cdata = cluster_clean(cdata, distance, min_neighboor)
xf = [cdata.red.xf;cdata.orange.xf];
yf = [cdata.red.yf;cdata.orange.yf];
zf = [cdata.red.zf;cdata.orange.zf];
colors = [cdata.red.xf*0+1;cdata.orange.xf*0+2];
clusters = dbscan([xf,yf,zf],distance,min_neighboor);
remove_id = clusters == -1;


cdata = remove_cdata(cdata,'red',remove_id(colors == 1));
cdata = remove_cdata(cdata,'orange',remove_id(colors == 2));

end