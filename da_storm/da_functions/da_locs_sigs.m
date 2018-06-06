function [xf_out,xf_cout, yf_out,yf_cout,sx_out, sx_cout, sy_out, sy_cout, N_out, N_cout,off_out, off_cout, fnout, llv_out, y, inloc,xin, yin] = da_locs_sig(i2, fnum, cents, angle)
% A living function that handles localization analysis


[m,n] = size(i2);
% load('z_cal.mat');
% zcurve = single([xsig, ysig, Ax, Ay, Bx, By, dx, dy, -gx, gy]);
maxi = 5000;
numthreads = 200;
xf_out = [];
yf_out = [];
N_out = [];
off_out = [];
xf_cout = [];
yf_cout = [];
N_cout = [];
off_cout = [];
llv_out = [];
xin = [];
yin = [];
fnout = [];
inloc = [];
sx_out = [];
sx_cout = [];
sy_out = [];
sy_cout = [];
y = zeros(n,1);
rads = single(deg2rad(angle));
if n > maxi
    rounds = floor(n/maxi);
    lefo = mod(n,maxi);
    
    for i = 1:rounds
        i3 = i2(:,(i-1)*maxi +1 :i*maxi);
        fna = fnum((i-1)*maxi +1 :i*maxi);
        cen = cents((i-1)*maxi +1 :i*maxi, :);
        [m,n] = size(i3);
        [xf, xc, yf, yc, Np, Nc,sx, sxc, sy, syc, off, offc, lv] = full_chain_loc(i3, numthreads, rads);
            
        
%         qind = Nc == -1;
%         offc = -offc;
%         yc = -yc;
%         offc(qind) = -1;
%         yc(qind) = -1;
%         
        ind = find( xc < 1 & yc <1 & xc > 0 & yc >0 & Nc > 0 & offc > 0);
        
        fnout = [fnout;fna(ind)];
        xf_out = [xf_out;xf(ind)+cen(ind,1)];
        yf_out = [yf_out;yf(ind)+cen(ind,2)];
        xin = [xin ; cen(ind, 1)];
        yin = [yin ; cen(ind, 2)];
        sx_out = [sx_out;sx(ind)];
        sy_out = [sy_out;sy(ind)];
        sx_cout = [sx_cout;sxc(ind)];
        sy_cout = [sy_cout;syc(ind)];
        y(ind + (i-1)*maxi) = 1;
        N_out = [N_out ; Np(ind)];
        off_out = [off_out; off(ind)];
        xf_cout = [xf_cout; xc(ind)];
        yf_cout = [yf_cout; yc(ind)];        
        N_cout = [N_cout; Nc(ind)];
        off_cout = [off_cout; offc(ind)];
        llv_out = [llv_out;lv(ind)];
        inloc = [inloc,i3(:,ind)];
        %         ajn_wait(mean(t), i,rounds);
        %         disp('Localizing');
    end
    i3 = i2(:,i*maxi+1:end);
    fna = fnum(i*maxi+1:end);
    cen = cents(i*maxi+1:end, :);
            [m,n] = size(i3);
        [xf, xc, yf, yc, Np, Nc,sx, sxc, sy, syc, off, offc, lv] = full_chain_loc(i3, numthreads, rads);
        
%         qind = Nc == -1;
%         offc = -offc;
%         yc = -yc;
%         offc(qind) = -1;
%         yc(qind) = -1;
        
        ind = find( xc < 1 & yc <1 & xc > 0 & yc >0 & Nc > 0 & offc > 0);
        
        fnout = [fnout;fna(ind)];
        xf_out = [xf_out;xf(ind)+cen(ind,1)];
        yf_out = [yf_out;yf(ind)+cen(ind,2)];
        xin = [xin ; cen(ind, 1)];
        yin = [yin ; cen(ind, 2)];
        sx_out = [sx_out;sx(ind)];
        sy_out = [sy_out;sy(ind)];
        sx_cout = [sx_cout;sxc(ind)];
        sy_cout = [sy_cout;syc(ind)];
        y(ind + (i-1)*maxi) = 1;
        N_out = [N_out ; Np(ind)];
        off_out = [off_out; off(ind)];
        xf_cout = [xf_cout; xc(ind)];
        yf_cout = [yf_cout; yc(ind)];
        N_cout = [N_cout; Nc(ind)];
        off_cout = [off_cout; offc(ind)];
        llv_out = [llv_out;lv(ind)];
        inloc = [inloc,i3(:,ind)];
        
else
    try
        fna = fnum;
        cen = cents;
            [m,n] = size(i2);
       [xf, xc, yf, yc, Np, Nc,sx, sxc, sy, syc, off, offc, lv] = full_chain_loc(i2, numthreads, rads);
                
        
        % correct minus sign errors for CRLB
%         qind = Nc == -1;
%         offc = -offc;
%         yc = -yc;
%         offc(qind) = -1;
%         yc(qind) = -1;
%         
    
        ind = find(xc < 1 & yc <1 & xc > 0 & yc >0 & Nc > 0 & offc > 0);
        
        fnout = [fnout;fna(ind)];
        xf_out = [xf_out;xf(ind)+cen(ind,1)];
        yf_out = [yf_out;yf(ind)+cen(ind,2)];
        xin = [xin ; cen(ind, 1)];
        yin = [yin ; cen(ind, 2)];
        sx_out = [sx_out;sx(ind)];
        sy_out = [sy_out;sy(ind)];
        sx_cout = [sx_cout;sxc(ind)];
        sy_cout = [sy_cout;syc(ind)];
        y(ind) = 1;
        N_out = [N_out ; Np(ind)];
        off_out = [off_out; off(ind)];
        xf_cout = [xf_cout; xc(ind)];
        yf_cout = [yf_cout; yc(ind)];
%         zf_cout = [zf_cout; zc(ind)];
        N_cout = [N_cout; Nc(ind)];
        off_cout = [off_cout; offc(ind)];
        llv_out = [llv_out;lv(ind)];
        inloc = [inloc,i2(:,ind)];
    catch lsterr
    end
end

disp('Done Localizing');
end


