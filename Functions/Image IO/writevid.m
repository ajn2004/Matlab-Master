function writevid(ims, fname, fps)
v = VideoWriter(fname, 'Motion JPEG AVI');
v.Quality = 98;
v.FrameRate = fps;
[m,n,o] = size(ims);
open(v);
for i = 1:o
    writeVideo(v,ims(:,:,i)./max(max(ims(:,:,i))));
end
close(v);