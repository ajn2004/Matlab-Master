function [xf_all, yf_all, N_all, sigx_all, sigy_all, off_all, framenum_all, number, xf_crlb, yf_crlb, N_crlb, sigx_crlb, sigy_crlb, off_crlb,llv, points] = func_time_series_spt(ave3, psize, sigma2)
close all
takeout = [];

ifin = ave3.*(ave3>0);

%% Manual Select of measuring regions
ishow = std(ifin,1,3);
[points] = select_cents(ishow,'jet', psize); % select points
count = 1;
xf_all = [];
yf_all = [];
framenum_all = [];
sigx_all = [];
sigy_all = [];
off_all = [];
N_all = [];
number = [];
xf_crlb = [];
yf_crlb = [];
N_crlb = [];
sigx_crlb = [];
sigy_crlb = [];
off_crlb = [];
number = [];
% section and fit
for i = 1:numel(points(:,1))
    % cut out regions for subsequent analysis
    [i2s, framenum, beta0] = cut_them_up(ave3, points(i,:), psize);
    for j = 1:numel(i2s(1,1,:))
        [fits, crlbs, llv] = func_mle_crlb(i2s(:,:,j), beta0(j,1), beta0(j,2), sigma2);
            xf = fits(1);
            yf = fits(2);
             N = fits(3);
          sigx = fits(4);
          sigy = fits(5);
        offset = fits(6);
        
        % check that everything makes sense
        if xf == xf && yf == yf && N == N && sigx == sigx && sigy == sigy && offset == offset % are the numbers real?
            if abs(xf) <= psize && abs(yf) <= psize && sigy > 0 && sigx > 0  % do the numbers make sense?
                if abs(xf) ~= inf && abs(yf) ~= inf && abs(N) ~= inf && abs(offset) ~= inf && abs(sigy) ~= inf && abs(sigx) ~= inf % make sure you aren't dealing with infinity
                    xf_all = [xf_all ; xf];
                    yf_all = [yf_all ; yf];
                    framenum_all = [framenum_all ; framenum(j)];
                    sigx_all = [sigx_all ; sigx];
                    sigy_all = [sigy_all ; sigy];
                    off_all = [off_all ; offset];
                    N_all = [N_all ; N];
                    number = [number ; i];
                    xf_crlb = [xf_crlb ; crlbs(1)];
                    yf_crlb = [yf_crlb ; crlbs(2)];
                    N_crlb = [N_crlb ; crlbs(3)];
                    sigx_crlb = [sigx_crlb ; crlbs(4)];
                    sigy_crlb = [sigy_crlb ; crlbs(5)];
                    off_crlb = [off_crlb ; crlbs(6)];
                end
            end
        end
    end
end
