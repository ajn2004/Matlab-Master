function func_batch_h2_tol_ps(fname)
% batch routine for files
try
load('C:\Users\ajnel\Documents\GitHub\Matlab-Master\Hurricane\Tolfile.mat');
catch
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\Tolfile.mat');
end
load(fname,'cdata','cal','pixw','q');

% Identify number of colors we're working with
flag = 0;
try
    xf = cdata.red.xf; % if red loads successfully, start red tolerances
    flag = 1; % if red exists assume 2-color, if orange also exists flag == 1
    
catch lsterr
end

try
    xf = cdata.orange.xf; % if orange loads successful, start orange tolerances
    if flag == 0 % if no red data this is an 'orange only'
        flag = 2;
    end
    
catch
    flag = 3; % Orange channel failed indicates must have red data, or will cause error
end
channel_flag = flag;
% channel_flag
% if channel flag = 1, dual color
% if channel flag = 2, orange only
% if channel flag = 3, red only

%% defining statistical red
if channel_flag == 1 || channel_flag == 3 % load red tolerances
    %Eliminate out of view molecue
    fnames = fieldnames(cdata.red);
    ind = cdata.red.xf < 0 | cdata.red.xf > 200 |cdata.red.yf < 0 | cdata.red.yf > 200;
    for i = 1:numel(fnames)
        cdata.red.(fnames{i})(ind,:) = [];
    end
    %Define fractionals
    cdata.red.fr_N =  cdata.red.crlbs(:,3).^0.5./cdata.red.fits(:,3);
    cdata.red.fr_sx = cdata.red.crlbs(:,4).^0.5./cdata.red.fits(:,4);
    cdata.red.fr_sy = cdata.red.crlbs(:,5).^0.5./cdata.red.fits(:,5);
    cdata.red.fr_o =  cdata.red.crlbs(:,6).^0.5./cdata.red.fits(:,6);
    cdata.red.ilv = cdata.red.llv(:)./cdata.red.fits(:,3);
    cdata.red.eps = abs(cdata.red.fits(:,4)./cdata.red.fits(:,5));
    cdata.red.snr =(cdata.red.fits(:,3)./(cdata.red.fits(:,3)+(2*6+1)^2*cdata.red.fits(:,6)).^0.5);
    %% Apply  red Tolerance
    
    % Build Red
    ind = cdata.red.fits(:,3) > tol.r.N_tol(1) & cdata.red.fits(:,3) < tol.r.N_tol(2); % Photon Tolerance
    ind = ind & cdata.red.snr >= tol.r.minsnr; % Signal to Noise tolerances
    ind = ind & abs(cdata.red.framenumber - mean(tol.r.flims)) <= diff(tol.r.flims)/2; % Framenumber tolerances
    ind = ind & abs(cdata.red.zf*q-mean(tol.r.zlims)) <= diff(tol.r.zlims)/2; % axial tolerances
    ind = ind & abs(cdata.red.fits(:,6)-mean(tol.r.offlim)) <= diff(tol.r.offlim)/2;   % offset tolerances
    ind = ind & (abs(cdata.red.fits(:,4)).*abs(cdata.red.fits(:,5))).^0.5 > tol.r.s_tol(1) & (abs(cdata.red.fits(:,4)).*abs(cdata.red.  fits(:,5))).^0.5 < tol.r.s_tol(2); % sigma Tolerance
    ind = ind & q*cdata.red.crlbs(:,1).^.5 < tol.r.lat_max & q*cdata.red.crlbs(:,2).^.5 < tol.r.lat_max; % Lateral Uncertainty Tolerance
    ind = ind & cdata.red.ilv > tol.r.iln; % llv tolerance
    ind = ind & cdata.red.fr_N < tol.r.frac_lim & abs(cdata.red.fr_o) < tol.r.off_frac; % Fraction photon tolerance
    ind = ind & cdata.red.fr_sx < tol.r.frac_lim & cdata.red.fr_sy < tol.r.frac_lim; % Fraction width tolerance
    rind = ind; % Index represents all molecules that PASS tolerances
    fnames = fieldnames(cdata.red);
    ind = logical(1-rind); % to remove failed localizations we want 1-pass
    for i = 1:numel(fnames) % remove data points
        cdata.red.(fnames{i})(ind,:) = [];
    end
end
%%  defining statistical orange
if channel_flag == 1 || channel_flag == 2 % load orange tolerances
    cdata.orange.xf = cdata.orange.xf(:);
    cdata.orange.yf = cdata.orange.yf(:);
    ind = cdata.orange.xf < 0 | cdata.orange.xf > 200 | cdata.orange.yf < 0 | cdata.orange.yf > 200;
    
    fnames = fieldnames(cdata.orange);
    
    for i = 1:numel(fnames)
        cdata.orange.(fnames{i})(ind,:) = [];
    end
    cdata.orange.fr_N =  cdata.orange.crlbs(:,3).^0.5./cdata.orange.fits(:,3);
    cdata.orange.fr_sx = cdata.orange.crlbs(:,4).^0.5./cdata.orange.fits(:,4);
    cdata.orange.fr_sy = cdata.orange.crlbs(:,5).^0.5./cdata.orange.fits(:,5);
    cdata.orange.fr_o =  cdata.orange.crlbs(:,6).^0.5./cdata.orange.fits(:,6);
    cdata.orange.ilv = cdata.orange.llv(:)./cdata.orange.fits(:,3);
    cdata.orange.eps = abs(cdata.orange.fits(:,4)./cdata.orange.fits(:,5));
    cdata.orange.snr =(cdata.orange.fits(:,3)./(cdata.orange.fits(:,3)+(2*6+1)^2*cdata.orange.fits(:,6)).^0.5);
    
    %% Apply  orange Tolerance
    % Build Orange
    ind = cdata.orange.fits(:,3) > tol.o.N_tol(1) & cdata.orange.fits(:,3) < tol.o.N_tol(2); % Photon Tolerance
    ind = ind & cdata.orange.snr >= tol.o.minsnr; % Signal to Noise tolerances
    ind = ind & abs(cdata.orange.framenumber - mean(tol.o.flims)) <= diff(tol.o.flims)/2; % Framenumber tolerances
    ind = ind & abs(cdata.orange.zf*q-mean(tol.o.zlims)) <= diff(tol.o.zlims)/2; % axial tolerances
    ind = ind & abs(cdata.orange.fits(:,6)-mean(tol.o.offlim)) <= diff(tol.o.offlim)/2;   % offset tolerances
    ind = ind & (abs(cdata.orange.fits(:,4)).*abs(cdata.orange.fits(:,5))).^0.5 > tol.o.s_tol(1) & (abs(cdata.orange.fits(:,4)).*abs(cdata.orange.  fits(:,5))).^0.5 < tol.o.s_tol(2); % sigma Tolerance
    ind = ind & q*cdata.orange.crlbs(:,1).^.5 < tol.o.lat_max & q*cdata.orange.crlbs(:,2).^.5 < tol.o.lat_max; % Lateral Uncertainty Tolerance
    ind = ind & cdata.orange.ilv > tol.o.iln; % llv tolerance
    ind = ind & cdata.orange.fr_N < tol.o.frac_lim & abs(cdata.orange.fr_o) < tol.o.off_frac; % Fraction photon tolerance
    ind = ind & cdata.orange.fr_sx < tol.o.frac_lim & cdata.orange.fr_sy < tol.o.frac_lim; % Fraction width tolerance
    oind = ind; % Index represents all molecules that PASS tolerances
    fnames = fieldnames(cdata.orange);
    ind = logical(1-oind);
    for i = 1:numel(fnames)
        cdata.orange.(fnames{i})(ind,:) = [];
    end
end

% Scream here for sanity
save([fname(1:end-4),'_tol.mat']);
end