function dps = das_peaks(ipca, varargin)
% Das Peaks returns an image DPS the same size as iprod that contains 1s at
% sites where a molecule likely exists

%% varargin stuff
if numel(varargin) < 1
    thresh = 20;
elseif numel(varargin) == 1
    thresh = varargin{1};
end
%% post-varargin
[m,n,p] = size(ipca);
dps = ipca.*0;
disp('Finding your Molecules');
for ii = 1:p
    im1 = ipca(:,:,ii);
    [sim1, I] = sort(im1(:));
    
    % Look at pixel values of the 50 - 75% percentiles of pixels
    x = round(numel(sim1)*0.5):round(numel(sim1)*0.75);
    y = sim1(x).';
    a = polyfit(x,y,1);  % create a linear fit to subtract off
    xp = round(numel(sim1)*0.5):numel(sim1); % create lin space for all pixels
    yp = a(1)*xp + a(2);
    di1 = sim1(xp)-yp.';
%     di1 = sim1;
    ind = xp(1) + find(di1 > thresh*mean(di1(end-5:end))/100, 1, 'first');
    Inds = I(ind:end);
    clear x y a sim1 xp yp di1 ind im1
    
    rows = mod(Inds,m) + 1;
    cols = floor(Inds/m) +1;

    z1 = zeros(m,n);
    for i = 1:numel(rows)
        z1(rows(i),cols(i)) = 1;
    end

    x = [];
    y = [];
    for j = 1:5
        s = strel('disk',6-j);
        z2 = imerode(z1,s);
        M = regionprops(logical(z2),'Centroid');
        % clear s z1 z2
        xs = [];
        ys =[];
        
        for i = 1:numel(M)
            xs = [xs;M(i).Centroid(1)];
            ys = [ys;M(i).Centroid(2)];
        end
        
        for i = 1:numel(x)
            dist = ((xs - x(i)).^2 + (ys-y(i)).^2).^0.5;
            ind = find (dist <= 3.5);
            xs(ind) = [];
            ys(ind) = [];
        end
        
        x = [x;xs];
        y = [y;ys];
    end
    for i = 1:numel(x)
        dps(round(y(i)),round(x(i)),ii) = 1;
    end
    clear x y xs ys z2 M s z1 rows cols clims
end
