function cfits = cen_locs(ims,fnum,cents)

[m,n,o] = size(ims);
pixw = (m-1)/2;
[X,Y] = meshgrid(-pixw:pixw);
for i = 1:o
    cfits(i,1) = sum(sum((X.*ims(:,:,i))/sum(sum(ims(:,:,i)))));
    cfits(i,2) = sum(sum((Y.*ims(:,:,i))/sum(sum(ims(:,:,i)))));
    cfits(i,3) = sum(sum(ims(:,:,i)));
end
cfits(:,1) = cfits(:,1) + cents(:,1);
cfits(:,2) = cfits(:,2) + cents(:,2);