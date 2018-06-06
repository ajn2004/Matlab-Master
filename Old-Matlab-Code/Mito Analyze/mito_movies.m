%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mito movies
%
% A workspace to make movies and figures for Jaime's lab meeting
% AJN 1/31/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;  % preclean line
roi = 5;
rec_movie = 'y'; % record movies?
% detect and find file
finfo = dir('*mito.mat');

% if no file detected ask for user assistance
if isempty(finfo)
    [fname, fpath] = uigetfile('*mito.mat');  % UI for selecting file
    load([fpath,fname]);  %load file into workspace
    cd(fpath); % change directory to chosen file's directory
else
    load(finfo.name);
end
% Everything we care about is wrapped up in this one variable 'fin'

mkdir('Analysis');  % make a directory to save stuff in

f1 = figure;
% Figure of red + green to demonstrate need for analysis
imy = fin(roi).green(:,:,1:3)*0;
imy(:,:,1) = fin(roi).red(:,:,1);
subplot(2,2,3);
imagesc(imy./(max(imy(:))));
axis image
title('Red Channel');
imy(:,:,2) = fin(roi).green(:,:,1);
subplot(2,2,[1,2])
imagesc(imy./(max(imy(:))))
axis image
title('Merge');
imy(:,:,1) = imy(:,:,1)*0;
subplot(2,2,4);
imagesc(imy./(max(imy(:))));
axis image
title('Green Channel');

% This figure will be used to generate the movie of mito ma
f2 = figure('Position',[244 154 915 717]);
[m,n,o] = size(fin(roi).green);
for i = 1:o % loop over frames
    subplot(2,3,1);
    imagesc(fin(roi).green(:,:,i));
    axis image
    title('Green');
    subplot(2,3,2);
    imagesc(fin(roi).mask(:,:,i));
    axis image
    title('Mask');
    subplot(2,3,3);
    imagesc(fin(roi).mask(:,:,i).*fin(roi).green(:,:,i));
    axis image
    title('Green Signal');
    
    subplot(2,3,4);
    imagesc(fin(roi).red(:,:,i));
    axis image
    title('Red');
    subplot(2,3,5);
    imagesc(fin(roi).mask(:,:,i));
    axis image
    title('Mask');
    subplot(2,3,6);
    imagesc(fin(roi).mask(:,:,i).*fin(roi).red(:,:,i));
    axis image
    title('Red signal');
%     colormap('gray')
    M(i) = getframe(gcf);
    
end
if strcmp(rec_movie,'y') || strcmp(rec_movie,'Y')
    movie2gif(M , [finfo(1).name(1:end-4),'_roi_',num2str(roi),'.gif'], 'DelayTime', 0.1);
end
ig = readtiff('c1_d1_10ms_2mM_1_btc_green.tif');
ir = readtiff('c1_d1_10ms_2mM_1_btc_red.tif');

igs = ig(:,:,10);

f3 = figure
imagesc(igs);
axis image;
title('Raw Mito Image')

rigs = rollingball(igs);

f4 = figure
imagesc(rigs);
axis image
title('Rolling Ball Subtraction Mito');

% f5 = figure
% irigs = ImsPCA(rigs,3);
% imagesc(irigs);
% axis image
% title('PCA Processed Mito');

f6 = figure
imagesc(rigs.*(rigs > 0.1*max(rigs(:))))
axis image
title('Threshold Mito');
xlabel('Threshold set to 10% max pixel value')


%% Mask Exclusion Figures
f7 = figure('Position',[560 344 560 604])
[m,n,o] = size(fin(roi).red);
for i = 1:o
    subplot(2,2,1);
    irs = fin(roi).red(:,:,i);
    msk = fin(roi).mask(:,:,i);
    msk2 = fin(roi).amask(:,:,i);
    imagesc(irs);
    axis image
    title('Raw Red Channel');
    subplot(2,2,2);
    imagesc(irs.*msk);
    axis image
    title('Mito Red Channel');
    subplot(2,2,3);
    imagesc(irs.*msk2);
    axis image
    title('Cell Red Channel');
    subplot(2,2,4);
    imagesc(irs.*logical(1-msk2));
    axis image
    title('Non-Cell Red Channel')
    drawnow
    L(i) = getframe(gcf);
end


%% Helper Functions
function movie2gif(mov, gifFile, varargin)
% ==================
% Matlab movie to GIF Converter.
%
% Syntax: movie2gif(mov, gifFile, prop, value, ...)
% =================================================
% The list of properties is the same like for the command 'imwrite' for the
% file format gif:
%
% 'BackgroundColor' - A scalar integer. This value specifies which index in
%                     the colormap should be treated as the transparent
%                     color for the image and is used for certain disposal
%                     methods in animated GIFs. If X is uint8 or logical,
%                     then indexing starts at 0. If X is double, then
%                     indexing starts at 1.
%
% 'Comment' - A string or cell array of strings containing a comment to be
%             added to the image. For a cell array of strings, a carriage
%             return is added after each row.
%
% 'DelayTime' - A scalar value between 0 and 655 inclusive, that specifies
%               the delay in seconds before displaying the next image.
%
% 'DisposalMethod' - One of the following strings, which sets the disposal
%                    method of an animated GIF: 'leaveInPlace',
%                    'restoreBG', 'restorePrevious', or 'doNotSpecify'.
%
% 'LoopCount' - A finite integer between 0 and 65535 or the value Inf (the
%               default) which specifies the number of times to repeat the
%               animation. By default, the animation loops continuously.
%               For a value of 0, the animation will be played once. For a
%               value of 1, the animation will be played twice, etc.
%
% 'TransparentColor' - A scalar integer. This value specifies which index
%                      in the colormap should be treated as the transparent
%                      color for the image. If X is uint8 or logical, then
%                      indexing starts at 0. If X is double, then indexing
%                      starts at 1
%
% *************************************************************************
% Copyright 2007-2013 by Nicolae Cindea.

if (nargin < 2)
    error('Too few input arguments');
end

if (nargin == 2)
    frameNb = size(mov, 2);
    isFirst = true;
    h = waitbar(0, 'Generate GIF file...');
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, ~] = frame2im(mov(i));
        if (exist('rgb2ind', 'file'))
            [IND, map] = rgb2ind(RGB,256);
        else
            [IND, map] = aRGB2IND(RGB);
        end
        if isFirst
            imwrite(IND, map, gifFile, 'gif');
            isFirst=false;
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append');
        end
    end
    close(h);
end

if (nargin > 2)
    h = waitbar(0, 'Generate GIF file...');
    frameNb = size(mov, 2);
    isFirst = true;
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, ~] = frame2im(mov(i));
        if (exist('rgb2ind', 'file'))
            [IND, map] = rgb2ind(RGB,256);
        else
            [IND, map] = aRGB2IND(RGB);
        end
        if isFirst
            args = varargin;
            imwrite(IND, map, gifFile, 'gif', args{:});
            isFirst=false;
            
            % supress 'LoopCount' option from the args!!
            args = varargin;
            l = length(args);
            
            posLoopCount = 0;
            for ii = 1:l
                if(ischar(args{ii}))
                    if strcmp(args{ii}, 'LoopCount')
                        posLoopCount = ii;
                    end
                end
            end
            if (posLoopCount)
                args = {args{1:posLoopCount-1}, args{posLoopCount+2:end}};
            end
            
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append', ...
                args{:});
        end
    end
    close(h);
end
