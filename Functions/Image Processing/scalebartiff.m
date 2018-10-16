% scalebar tiff
% add a scale bar to a tiff series
mkdir('Bar');
scb = 5; % um
q = 0.133; % um/pix
sxbp = round(scb/q); % scb in pix
sch = 5; % pix
xst = 3;% pix
yst = 3;% pix

% [fname, fpath] = uigetfile('*tif');
files = dir('*tif');
for k = 1:numel(files)
    im1 = readtiff(files(k).name);
    [m,n,o] = size(i1);
    i1 = ip1;
    for i = 1:o
        maxi = max(max(i1(:,:,i)));
        mini = min(min(i1(:,:,i)));
        maxi = (maxi + mini)/2;
        i1(yst:yst + sch, xst:xst + sxbp, i) = maxi;
        %     imagesc(im1(:,:,i));
        %     drawnow
    end
    imagesc(i1(:,:,10));
    writetiff(mean(im1,3),['Bar\',files(k).name(1:end-4),'_bar']);
end