function [im1] = read_im(str)
im1 = [];
% fname = [];
% fpath = [];
ext = str(end-3:end);
% disp(ext)

if strcmp(ext,'.tif')
    im1 = readtiff(str);
elseif strcmp(ext,'fits')
    im1 = fitsread(str);
else
    disp('Extension not recognized');
end


% switch str
%     case '.tif'
% %         [im1, fname, fpath] = readtiff();
%         im1 = readtiff(str);
%     case 'fits'
% %         [fname, fpath] = uigetfile('*fits');
% %         im1 = fitsread([fpath,fname]);
%         im1 = fitsread(str);
%     otherwise
%         disp('Extension not recognized');
% end