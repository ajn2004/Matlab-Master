% make it a gif
function make_it_a_gif(i1,fname, varargin)
if nargin == 3
    clims = varargin{1};
end
[m,n,o] = size(i1);
cmap = colormap(jet(200));

for i = 1:o
    if nargin < 3
        im = color_it(i1(:,:,i),cmap);
    else
        im = color_it(i1(:,:,i),cmap,clims);
    end
    [imind, cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,[fname,'.gif'],'gif','LoopCount',Inf,'DelayTime',0.03);
    else
        imwrite(imind,cm,[fname,'.gif'],'gif', 'WriteMode', 'append','DelayTime',0.03);
    end
end