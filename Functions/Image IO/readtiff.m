function [im2, fname, fpath] = readtiff(varargin)
%READTIFF Read a tiff image faster than imread
%   A = READTIFF(FILENAME) reads a tiff from the file specified by the
%   string FILENAME. FILENAME must be in the current directory, or have the
%   path specified. A will be a matrix of size image width x image height x
%   frames in image. Without specifiying N_END and N_START all frames in
%   the files will be read.
%
%   A = READTIFF(FILENAME,N_START, N_END) reads a tiff from the file specified by the
%   string FILENAME. FILENAME must be in the current directory, or have the
%   path specified. A will be a matrix of size image width x image height x
%   N_END - N_START + 1. Specifying N_START and N_END will load only frames
%   starting at N_START and ending at N_END. Setting N_END to 0 will load
%   all frames from N_START to the last frame
im1 = [];
% find number of arguments after fname
numvar = length(varargin);
w = warning ('off','all');
if numvar == 0 
    [fname, fpath] = uigetfile('*.tif');
else
    fname = varargin{1};
end
if numvar > 3
    error('Too many optional inputs, please only specify n_start and n_end')
elseif numvar < 2
%     disp('Not enough specification variables, starting from 1 to finish');
    n_start = 1;
    n_end = 0;
elseif numvar == 2
    n_start = varargin{2};
    n_end = varargin{2};
elseif numvar == 3
    if varargin{2} > varargin{3}
        n_start = varargin{3};
        n_end = varargin{2};
    else
        n_start = varargin{2};
        n_end = varargin{3};
    end
end

InfoImage=imfinfo(fname);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
if n_end == 0
    n_end = NumberImages;
end
im1=zeros(nImage,mImage,n_end - n_start +1,'uint16');
TifLink = Tiff(fname, 'r');
count = 1;
for i=n_start:n_end
%     tic
   TifLink.setDirectory(i);
   im1(:,:,count)=TifLink.read();
%    t(count) = toc;
   count= count +1;   
%    ajn_wait(mean(t), i, (n_end-n_start +1));
%    disp('Loading');
end
TifLink.close();
w = warning ('on','all');
im2 = single(im1);
end
