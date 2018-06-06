% scalebar tiff
% add a scale bar to a tiff series
mkdir('Bar');
scb = 5; % um
q = 0.133; % um/pix
sxbp = round(scb/q); % scb in pix
sch = 5; % pix
xst = 10;% pix
yst = 10;% pix

% [fname, fpath] = uigetfile('*tif');
files = dir('*tif');
for k = 1:numel(files)
    im1 = readtiff(files(k).name);
    [m,n,o] = size(im1);
    
    for i = 1:o
        maxi = max(max(im1(:,:,i)));
        im1(yst:yst + sch, xst:xst + sxbp, i) = maxi;
        %     imagesc(im1(:,:,i));
        %     drawnow
    end
    
    writetiff(mean(im1,3),['Bar\',files(k).name(1:end-4),'_bar']);
end