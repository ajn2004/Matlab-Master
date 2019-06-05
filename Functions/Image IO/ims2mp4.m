function ims2mp4(i1,fname,varargin)
% This function will convert a image stack to an MP4 sequence
switch nargin
    case 0
        error('Need 2 inputs (images, name)');
    case 1
        error('Need 2 inputs (images, name)');
    case 2
        fps = 30;
        qual = 75;
    case 3
        fps = varargin{1};
        qual = 75;
    case 4
        fps = varargin{1};
        qual = varargin{2};
end
if isempty(strfind(fname,'avi'))
    fname = [fname,'.avi'];
end
v = VideoWriter(fname, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = qual;
[m,n,o] = size(i1);
open(v)
for i = 1:o
    tic
    
    writeVideo(v,i1(:,:,i));
    t(i) = toc;
    ajn_wait(t,i,o)
end
close(v)
