
3
.% scalebar tiff
% add a scale bar to a tiff series
mkdir('Bar');
scb = 10; % um
q = 0.128; % um/pix
sxbp = round(scb/q); % scb in pix
sch = 5; % pix
xst = 30;% pix
yst = 3;% pix

% [fname, fpath] = uigetfile('*tif');
files = dir('*tif');
% for k = 1:numel(files)
%     im1 = readtiff(files(k).name);
    [m,n,o] = size(im1);
%     im1 = ip1;
    for i = 1:o
        maxi = max(max(im1(:,:,i)));
        mini = min(min(im1(:,:,i)));
        maxi = (maxi + mini)/2;
        im1(yst:yst + sch, xst:xst + sxbp, i) = maxi;
        %     imagesc(im1(:,:,i));
        %     drawnow
    end
%     imagesc(i1(:,:,10));
%     writetiff(mean(im1,3),['Bar\',files(k).name(1:end-4),'_bar']);
% end