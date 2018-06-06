clearvars
% close all
%

i1 = double(readtiff('storm.tif'));
[m,n,o] = size(i1);
pad = 1;
iprod = rollingball(i1);
count = 1;
li = randperm(o);
% for k = 1:10
    i2 = [zeros(pad,n+2*pad);zeros(m,pad),iprod(:,:,200),zeros(m,pad);zeros(pad,n+2*pad)];
    
    tic
    % image filter
    for i = pad +1:m+pad
        for j = pad+1:n+pad
            i3 = i2(i-pad:i+pad,j-pad:j+pad);
            ims(count,:) = i3(:).';
            count = count +1;
        end
    end
% end
toc
for i = 1:numel(ims(1,:))
    ims(:,i) = ims(:,i) - mean(ims(:,i));
end
close all
thre = 0.5;
[V, D] = eig(cov(ims));
plot(max(D))
dthresh = thre*max(D(:));
[row, col] = find(D>dthresh);
k = 1
Vproj = V(:,col);
vals = Vproj.'*ims.';
% while true
for i = 1:m
    for j= 1:n
        i4(i,j) = vals((i-1)*n +j);
    end
end

imagesc(i4)
figure
surf(i4)
title('i4')
figure
surf(i1);
title('i1')