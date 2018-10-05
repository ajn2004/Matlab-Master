function [xfo, xfc, yfo, yfc, N, Nco, Zo, Zco, fnout] = slim_locs(i2,fnum, cents, angle, zin, zcin, nin, ncin)
% Slim Locs is a script to analyze data w/ the gaussian fit after
% performing the cspline fit to get Z data. 
maxi = 5000;
[m,n] = size(i2);
i2 = double(i2);
rads = double(angle);
xfo = [];
xfc = [];
yfo = [];
yfc = [];
N = [];
Nco = [];
Zo = [];
Zco = [];
fnout = [];
if n > maxi
    rounds = floor(n/maxi);
    lefo = mod(n,maxi);
    
    for i = 1:rounds
        i3 = i2(:,(i-1)*maxi +1 :i*maxi);
        fna = fnum((i-1)*maxi +1 :i*maxi);
        cen = cents((i-1)*maxi +1 :i*maxi, :);
        zs = zin((i-1)*maxi +1 :i*maxi);
        zcs = zcin((i-1)*maxi +1 :i*maxi);
        ns = nin((i-1)*maxi +1 :i*maxi);
        ncs = ncin((i-1)*maxi +1 :i*maxi);
        [m,n] = size(i3);
        switch m % Switch localization algorithm based on size
            case 9
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_3(i3, 100, rads, 15);
            case 25
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_5(i3, 100, rads, 15);
            case 49
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_7(i3, 100, rads, 15);
            case 81
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_9(i3, 100, rads, 15);
            case 121
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_11(i3, 100, rads, 15);
            case 169
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_13(i3, 100, rads, 15);
            case 225
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_15(i3, 100, rads, 15);
                
            otherwise
        end
        
        yc = -yc;
        ind = Nc ~= -1;
        fnout=[fnout;fna(ind)];
        xfo = [xfo;xf(ind)+cen(ind,1)];
        yfo = [yfo;yf(ind)+cen(ind,2)];
        xfc = [xfc; xc(ind)];
        yfc = [yfc; yc(ind)];
        Zo = [Zo; zs(ind)];
        Zco = [Zco; zcs(ind)];
        N = [N;ns(ind)];
        Nco = [Nco; ncs(ind)];
        
        
    end
    i3 = i2(:,i*maxi+1:end);
    fna = fnum(i*maxi+1:end);
    cen = cents(i*maxi+1:end, :);
    
    zs = zin(i*maxi+1:end);
    zcs = zcin(i*maxi+1:end);
    ns = nin(i*maxi+1:end);
    ncs = ncin(i*maxi+1:end);
    [m,n] = size(i3);
    switch m % Switch localization algorithm based on size
        case 9
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_3(i3, 100, rads, 15);
        case 25
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_5(i3, 100, rads, 15);
        case 49
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_7(i3, 100, rads, 15);
        case 81
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_9(i3, 100, rads, 15);
        case 121
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_11(i3, 100, rads, 15);
        case 169
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_13(i3, 100, rads, 15);
        case 225
            [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_15(i3, 100, rads, 15);
            
        otherwise
    end
    
    yc = -yc;
    ind = Nc ~= -1;
    fnout=[fnout;fna(ind)];
    xfo = [xfo;xf(ind)+cen(ind,1)];
    yfo = [yfo;yf(ind)+cen(ind,2)];
    xfc = [xfc; xc(ind)];
    yfc = [yfc; yc(ind)];
    Zo = [Zo; zs(ind)];
    Zco = [Zco; zcs(ind)];
    N = [N;ns(ind)];
    Nco = [Nco; ncs(ind)];
    
else
    try
        fna = fnum;
        cen = cents;
        zs = zin;
        zcs = zcin;
        ns = nin;
        ncs = ncin;
        [m,n] = size(i2);
        switch m % Switch localization algorithm based on size
            case 9
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_3(i2, 100, rads, 15);
            case 25
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_5(i2, 100, rads, 15);
            case 49
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_7(i2, 100, rads, 15);
            case 81
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_9(i2, 100, rads, 15);
            case 121
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_11(i2, 100, rads, 15);
            case 169
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_13(i2, 100, rads, 15);
            case 225
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_15(i2, 100, rads, 15);
                
            otherwise
        end
        
        yc = -yc;
        ind = Nc ~= -1;
        fnout=[fnout;fna(ind)];
        xfo = [xfo;xf(ind)+cen(ind,1)];
        yfo = [yfo;yf(ind)+cen(ind,2)];
        xfc = [xfc; xc(ind)];
        yfc = [yfc; yc(ind)];
        Zo = [Zo; zs(ind)];
        Zco = [Zco; zcs(ind)];
        N = [N;ns(ind)];
        Nco = [Nco; ncs(ind)];
    catch lsterr
    end
end

end

