function [xf_out,xf_cout, yf_out,yf_cout, N_out, N_cout, sigx_out, sigx_cout,sigy_out, sigy_cout,off_out, off_cout, llv_out] = break_loc(im1)
[m,n,o] = size(im1);

[xg, yg] = (meshgrid(1:n,1:m));
xg = single(xg);
yg = single(yg);
i1 = filter_cut(im1(:,:,1),3).';
maxi = 20000;

xf_all = [];
yf_all = [];
N = [];
sigx_all = [];
sigy_all = [];
off_all = [];
xf_crlb = [];
yf_crlb = [];
N_crlb = [];
sigx_crlb = [];
sigy_crlb = [];
off_crlb = [];
llv = [];
xin = [];
yin = [];
if n*m > maxi
    rounds = floor(n*m/maxi);
    lefo = mod(m*n,maxi);
    
    for i = 1:rounds
        i2 = i1(:,(i-1)*maxi +1 :i*maxi);
        
        [xf,xc, yf,yc, Np,  Nc, sigx, sigxc, sigy, sigyc,off, offc, lv] = full_chain_loc(i2,100);
        
        ind = find(sigxc > 0 & xc < 1 & yc <1 & xc > 0 & yc >0);
        xf_all = [xf_all;xf(ind)+xg(ind)];
        yf_all = [yf_all;yf(ind)+yg(ind)];
        xin = [xin ; xg(ind)];
        yin = [yin ; yg(ind)];
        %           xf_all = [xf_all;xf(ind)];
        %           yf_all = [yf_all;yf(ind)];
        N = [N ; Np(ind)];
        sigx_all = [sigx_all;sigx(ind)];
        sigy_all = [sigy_all; sigy(ind)];
        off_all = [off_all; off(ind)];
        xf_crlb = [xf_crlb; xc(ind)];
        yf_crlb = [yf_crlb; yc(ind)];
        N_crlb = [N_crlb; Nc(ind)];
        sigx_crlb = [sigx_crlb; sigxc(ind)];
        sigy_crlb = [sigy_crlb; sigyc(ind)];
        off_crlb = [off_crlb; offc(ind)];
        llv = [llv;lv(ind)];
    end
    i2 = i1(:,i*maxi+1:end);
    [xf,xc, yf,yc, Np,  Nc, sigx, sigxc, sigy, sigyc,off, offc, lv] = full_chain_loc(i2,100);
    ind = find(sigxc > 0 & xc < 1 & yc <1 & xc > 0 & yc >0);
    xf_all = [xf_all;xf(ind)+xg(ind)];
    yf_all = [yf_all;yf(ind)+yg(ind)];
    xin = [xin ; xg(ind)];
    yin = [yin ; yg(ind)];
    %           xf_all = [xf_all;xf(ind)];
    %           yf_all = [yf_all;yf(ind)];
    N = [N ; Np(ind)];
    sigx_all = [sigx_all;sigx(ind)];
    sigy_all = [sigy_all; sigy(ind)];
    off_all = [off_all; off(ind)];
    xf_crlb = [xf_crlb; xc(ind)];
    yf_crlb = [yf_crlb; yc(ind)];
    N_crlb = [N_crlb; Nc(ind)];
    sigx_crlb = [sigx_crlb; sigxc(ind)];
    sigy_crlb = [sigy_crlb; sigyc(ind)];
    off_crlb = [off_crlb; offc(ind)];
    llv = [llv;lv(ind)];
else
    [xf,xc, yf,yc, Np,  Nc, sigx, sigxc, sigy, sigyc,off, offc, lv] = full_chain_loc(i1,100);
    ind = find(sigxc > 0 & xc < 1 & yc <1 & xc > 0 & yc >0);
    xf_all = [xf_all;xf(ind)+xg(ind)];
    yf_all = [yf_all;yf(ind)+yg(ind)];
    xin = [xin ; xg(ind)];
    yin = [yin ; yg(ind)];
    %           xf_all = [xf_all;xf(ind)];
    %           yf_all = [yf_all;yf(ind)];
    N = [N ; Np(ind)];
    sigx_all = [sigx_all;sigx(ind)];
    sigy_all = [sigy_all; sigy(ind)];
    off_all = [off_all; off(ind)];
    xf_crlb = [xf_crlb; xc(ind)];
    yf_crlb = [yf_crlb; yc(ind)];
    N_crlb = [N_crlb; Nc(ind)];
    sigx_crlb = [sigx_crlb; sigxc(ind)];
    sigy_crlb = [sigy_crlb; sigyc(ind)];
    off_crlb = [off_crlb; offc(ind)];
    llv = [llv;lv(ind)];
end

% Remove infs and NaNs
varys = {'yin','xin','xf_all','yf_all','N','sigx_all','sigy_all', 'off_all','off_crlb','xf_crlb','yf_crlb','N_crlb','sigy_crlb','sigx_crlb','llv'};
for i = 1:numel(varys)
    eval(['index = find(isnan(',varys{i},') | isinf(',varys{i},'));'])
    for j = 1:numel(varys)
        eval([varys{j},'(index) = [];']);
    end
end

% Remove Duplicates
% This is based on the idea that pixels are analyzed sequentially through
% column space, thus by finding overlapping localizations and selecting the
% best fit by the lowest combined fractional uncertainty
l = 1;
flags = [];
for i = 1:numel(xf_all)
    if ismember(i,flags)
    else
        ind = find( yin >= yin(i) -3 & yin <= yin(i) +3 &   xin >= xin(i) -3 & xin <= xin(i) + 3); % find entries within 2 pixels of localization
        frn = N_crlb(ind).^0.5./N(ind);
        fro = off_crlb(ind).^0.5./off_all(ind);
        frsx = sigx_crlb(ind).^0.5./sigx_all(ind);
        frsy = sigy_crlb(ind).^0.5./sigy_all(ind);
        bigf = frn + fro + frsx + frsy; %combine uncertainties
        indo = find(bigf == min(bigf)); %select minimum uncertainty
        
        xf_out(l,1) = xf_all(ind(indo(1)));
        xf_cout(l,1)= xf_crlb(ind(indo(1)));
        yf_out(l,1) = yf_all(ind(indo(1)));
        yf_cout(l,1) = yf_crlb(ind(indo(1)));
        N_out(l,1) = N(ind(indo(1)));
        N_cout(l,1) = N_crlb(ind(indo(1)));
        sigx_out(l,1) = sigx_all(ind(indo(1)));
        sigx_cout(l,1) = sigx_crlb(ind(indo(1)));
        sigy_out(l,1) = sigy_all(ind(indo(1)));
        sigy_cout(l,1) = sigy_crlb(ind(indo(1)));
        off_out(l,1) = off_all(ind(indo(1)));
        off_cout(l,1) = off_crlb(ind(indo(1)));
        llv_out(l,1) = llv(ind(indo(1)));
        l =l+1;
        ind(indo) = [];
        flags = [flags;ind];
    end
    
end
end

