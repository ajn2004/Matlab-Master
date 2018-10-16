function [xf_out,xf_cout, yf_out,yf_cout,zf_out,zf_cout, N_out, N_cout,off_out, off_cout, fnout, llv_out, its_out, cex_out, cey_out] = da_splines(iloc, fnum, cents, cal, pixw)

% Assume iloc from previous sections 
[m,n,o] = size(iloc); % get size of iloc
maxi = 1000; % set max images / gpu run
% ilocs = reshape(iloc,m^0.5,m^0.5,n); % reshape images
ilocs = iloc;
clear iloc
% Preallocate
xf_out = [];
yf_out = [];
zf_out = [];
N_out = [];
off_out = [];
xf_cout = [];
yf_cout = [];
zf_cout = [];
N_cout = [];
off_cout = [];
llv_out = [];
fnout = [];
its_out = [];
cex_out = [];
cey_out = [];


if o > maxi % If we can't localize it all at once, chunk it up
    rounds = floor(o/maxi);
    lefo = mod(o,maxi);
    
    for i = 1:rounds % Loop over chunks
        fna = fnum((i-1)*maxi +1 :i*maxi); % separate framenums
        cen = cents((i-1)*maxi +1 :i*maxi, :); % separate coms
        [P, CRLB, LogL] = mleFit_LM(ilocs(:,:,(i-1)*maxi +1 :i*maxi),5,50,single(cal.cspline.coeff), 0); % perform fit
%         ind = abs(P(:,1) - pixw ) <= 0.5 & abs(P(:,2) - pixw) <= 0.5;
        frx = CRLB(:,1).^0.5./P(:,1);
        fry = CRLB(:,2).^0.5./P(:,2);
%         ind = frx <= mean(frx) & fry <= mean(fry);
%         ind = P(:,1) > 6.5 & P(:,1) < 7 & P(:,2) > 3.25 & P(:,2) < 10.75;
    ind = 1:size(P(:,1));
        % Make final assignments
        fnout = [fnout;fna(ind)];
        xf_out = [xf_out;P((ind),1)+cen((ind),1)];
        yf_out = [yf_out;P((ind),2)+cen((ind),2)];
        N_out = [N_out ; P((ind),3)];
        off_out = [off_out; P((ind),4)];
        zf_out = [zf_out; (P((ind),5)-cal.cspline.z0)*cal.cspline.dz/1000]; % This puts zf_out in terms of um, we can convert to pixels on the outside
        its_out = [its_out; P((ind),6)];
        xf_cout = [xf_cout; CRLB((ind),1)];
        yf_cout = [yf_cout; CRLB((ind),2)];        
        N_cout = [N_cout; CRLB((ind),3)];
        off_cout = [off_cout; CRLB((ind),4)];
        zf_cout = [zf_cout; CRLB((ind),5)*(cal.cspline.dz/1000)^2]; % this converts the CRLB to um^2, we will convert to pixels on the outside
        llv_out = [llv_out;LogL(ind)];
    end
    
    % Finish the stack with a Try
    try
        
        fna = fnum(i*maxi+1:end); % separate framenums
        cen = cents(i*maxi+1:end, :); % separate coms
        [P, CRLB, LogL] = mleFit_LM(ilocs(:,:,i*maxi +1 :end),5,50,single(cal.cspline.coeff), 1); % perform fit
%         ind = abs(P(:,1) - pixw ) <= 0.5 & abs(P(:,2) - pixw) <= 0.5;
        frx = CRLB(:,1).^0.5./P(:,1);
        fry = CRLB(:,2).^0.5./P(:,2);
%         ind = frx <= mean(frx) & fry <= mean(fry);
        ind = 1:size(P(:,1));
        % Make final assignments
        fnout = [fnout;fna(ind)];
        xf_out = [xf_out;P((ind),1)+cen((ind),1)];
        yf_out = [yf_out;P((ind),2)+cen((ind),2)];
        N_out = [N_out ; P((ind),3)];
        off_out = [off_out; P((ind),4)];
        zf_out = [zf_out; (P((ind),5)-cal.cspline.z0)*cal.cspline.dz/1000]; % This puts zf_out in terms of um, we can convert to pixels on the outside
        its_out = [its_out; P((ind),6)];
        xf_cout = [xf_cout; CRLB((ind),1)];
        yf_cout = [yf_cout; CRLB((ind),2)];        
        N_cout = [N_cout; CRLB((ind),3)];
        off_cout = [off_cout; CRLB((ind),4)];
        zf_cout = [zf_cout; CRLB((ind),5)*(cal.cspline.dz/1000)^2]; % this converts the CRLB to um^2, we will convert to pixels on the outside
        llv_out = [llv_out;LogL(ind)];
    catch lsterr
    end
      
    
else
    try
        fna = fnum;
        cen = cents;
        [P, CRLB, LogL] = mleFit_LM(ilocs,5,50,single(cal.cspline.coeff), 1); % perform fit
        
%       ind = abs(P(:,1) - pixw ) <= 0.5 & abs(P(:,2) - pixw) <= 0.5; 
        frx = CRLB(:,1).^0.5./P(:,1);
        fry = CRLB(:,2).^0.5./P(:,2);
%         ind = frx <= mean(frx) & fry <= mean(fry);
%         ind = P(:,6) < 50;
      ind = 1:size(P(:,1));  % uncomment this line to catch all localizations
        % Make final assignments
        fnout = [fnout;fna(ind)];
        xf_out = [xf_out;P((ind),1)+cen((ind),1)];
        yf_out = [yf_out;P((ind),2)+cen((ind),2)];
        cex_out = [cex_out;cen((ind),1)];
        cey_out = [cey_out;cen((ind),2)];
        N_out = [N_out ; P((ind),3)];
        off_out = [off_out; P((ind),4)];
        zf_out = [zf_out; (P((ind),5)-cal.cspline.z0)*cal.cspline.dz/1000]; % This puts zf_out in terms of um, we can convert to pixels on the outside
        its_out = [its_out; P((ind),6)];
        xf_cout = [xf_cout; CRLB((ind),1)];
        yf_cout = [yf_cout; CRLB((ind),2)];        
        N_cout = [N_cout; CRLB((ind),3)];
        off_cout = [off_cout; CRLB((ind),4)];
        zf_cout = [zf_cout; CRLB((ind),5)*(cal.cspline.dz/1000)^2]; % this converts the CRLB to um^2, we will convert to pixels on the outside
        llv_out = [llv_out;LogL(ind)];
    catch lsterr
    end
end
    try
%         dx = floor(m/2);
%     xf_out = xf_out - dx;
%     yf_out = yf_out - dx;
    catch lsterr
    end
    disp('Done Localizing');
end
    
