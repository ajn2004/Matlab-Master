function [fits, crlbs, llv, fnout] = slim_locs(i2, varargin)
% Slim Locs is a script to analyze data w/ the gaussian fit after
% performing the cspline fit to get Z data.
% Variable input is expect to be of the form (i2, fnum, cents, rads,
% lpcount, thrds)

vars = nargin - 1;
[m,n,o] = size(i2);
for i = 1:o
    mini(i) = min(min(i2(:,:,i)));
    i2(:,:,i) = i2(:,:,i) - mini(i);
end
i2 = reshape(i2,m*n,o);
thrds = 100;
lpcnt = 10;
rads = 0;
if vars ==0
    fnum = 1:o; % supplement framenumber 
    cents = zeros(o,2);
end
if vars >= 1
    fnum = varargin{1};
end
if vars >= 2
    cents = varargin{2};
end
if vars >= 3
    rads = double(varargin{3});
end
if vars >= 4
     lpcnt = double(varargin{4});
end 
if vars >= 5
     thrds = double(varargin{5});
end 

maxi = 5000;
i2 = double(i2); % ensure double type variables
% Preallocate big variables
xfo = [];
xfc = [];
yfo = [];
yfc = [];
N = [];
Nco = [];
sxo = [];
sxco = [];
syo = [];
syco= [];
offo = [];
offco = [];
llv = [];
fnout = [];
n = o;
if n > maxi
    rounds = floor(n/maxi);    
    for i = 1:rounds
        i3 = i2(:,(i-1)*maxi +1 :i*maxi);
        fna = fnum((i-1)*maxi +1 :i*maxi);
        cen = cents((i-1)*maxi +1 :i*maxi, :);

        [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc(i3, thrds, rads, lpcnt);

         yc = -yc;
        ind = Nc ~= -1; % remove all failed localizations
      fnout = [fnout;fna(ind)];
        xs = cos(-rads)*xf - sin(-rads)*yf+cen(:,1);
      ys = sin(-rads)*xf + cos(-rads)*yf+cen(:,2);
        xfo = [xfo; xs(ind)]; % X-Y values will be rotated by rads, perform inverse transformation
        yfo = [yfo; ys(ind)];
        xfc = [xfc; xc(ind)];
        yfc = [yfc; yc(ind)];
          N = [N;Np(ind)];
        Nco = [Nco; Nc(ind)];
        sxo = [sxo;sx(ind)];
       sxco = [sxco;sxc(ind)];
        syo = [syo;sy(ind)];
       syco = [syco;syc(ind)];
       offo = [offo;off(ind)+mini(ind)];
      offco = [offco; offc(ind)];
        llv = [llv; lv(ind)];
        
    end
    i3 = i2(:,i*maxi+1:end);
    fna = fnum(i*maxi+1:end);
    cen = cents(i*maxi+1:end, :);
    if ~isempty(i3)
    [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc(i3, thrds, rads, lpcnt);
    
         yc = -yc;
        ind = Nc ~= -1; % remove all failed localizations
      fnout = [fnout;fna(ind)];
        xs = cos(-rads)*xf - sin(-rads)*yf+cen(:,1);
      ys = sin(-rads)*xf + cos(-rads)*yf+cen(:,2);
        xfo = [xfo; xs(ind)]; % X-Y values will be rotated by rads, perform inverse transformation
        yfo = [yfo; ys(ind)];
        xfc = [xfc; xc(ind)];
        yfc = [yfc; yc(ind)];
          N = [N;Np(ind)];
        Nco = [Nco; Nc(ind)];
        sxo = [sxo;sx(ind)];
       sxco = [sxco;sxc(ind)];
        syo = [syo;sy(ind)];
       syco = [syco;syc(ind)];
       offo = [offo;off(ind)];
      offco = [offco; offc(ind)];
        llv = [llv; lv(ind)];
    end
else
    try
        fna = fnum;
        cen = cents;
        
        [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc(i2, thrds, rads, lpcnt);
         yc = -yc;
        ind = Nc ~= -1; % remove all failed localizations
      fnout = [fnout;fna(ind)];
      xs = cos(-rads)*xf - sin(-rads)*yf+cen(:,1);
      ys = sin(-rads)*xf + cos(-rads)*yf+cen(:,2);
        xfo = [xfo; xs(ind)]; % X-Y values will be rotated by rads, perform inverse transformation
        yfo = [yfo; ys(ind)];
        xfc = [xfc; xc(ind)];
        yfc = [yfc; yc(ind)];
          N = [N;Np(ind)];
        Nco = [Nco; Nc(ind)];
        sxo = [sxo;sx(ind)];
       sxco = [sxco;sxc(ind)];
        syo = [syo;sy(ind)];
       syco = [syco;syc(ind)];
       offo = [offo;off(ind)];
      offco = [offco; offc(ind)];
        llv = [llv; lv(ind)];
    catch lsterr
    end
end
fits = [xfo(:,1),yfo(:,1),N,sxo,syo,offo]; %consolidate fits
crlbs = abs([xfc,yfc,Nco,sxco,syco,offco]); % consolidate errors
end

