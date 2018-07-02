function [im1, fname, fpath] = read_im(str)
im1 = [];
fname = [];
fpath = [];
switch str
    case 'tif'
        [im1, fname, fpath] = readtiff();
    case 'fits'
        [fname, fpath] = uigetfile('*fits');
        im1 = fitsread([fpath,fname]);
    otherwise
        disp('Extension not recognized');
end