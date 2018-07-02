function im1 = read_im(str)
im1 = [];
switch str
    case 'tif'
        im1 = readtiff();
    case 'fits'
        [fname, fpath] = uigetfile('*fits');
        im1 = fitsread([fpath,fname]);
    otherwise
        disp('Extension not recognized');
end