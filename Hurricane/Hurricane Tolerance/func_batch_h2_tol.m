function func_batch_h2_tol(fname)
% batch routine for files
load('Tolfile.mat');
load(fname,'cdata','cal','pixw','q');

%Define fractionals
cdata.red.fr_N =  cdata.red.crlbs(:,3).^0.5./cdata.red.fits(:,3);
cdata.red.fr_sx = cdata.red.crlbs(:,4).^0.5./cdata.red.fits(:,4);
cdata.red.fr_sy = cdata.red.crlbs(:,5).^0.5./cdata.red.fits(:,5);
cdata.red.fr_o =  cdata.red.crlbs(:,6).^0.5./cdata.red.fits(:,6);
cdata.red.ilv = cdata.red.llv(:)./cdata.red.fits(:,3);
cdata.red.eps = abs(cdata.red.fits(:,4)./cdata.red.fits(:,5));
cdata.red.snr =(cdata.red.fits(:,3)./(cdata.red.fits(:,3)+(2*6+1)^2*cdata.red.fits(:,6)).^0.5);

cdata.orange.fr_N =  cdata.orange.crlbs(:,3).^0.5./cdata.orange.fits(:,3);
cdata.orange.fr_sx = cdata.orange.crlbs(:,4).^0.5./cdata.orange.fits(:,4);
cdata.orange.fr_sy = cdata.orange.crlbs(:,5).^0.5./cdata.orange.fits(:,5);
cdata.orange.fr_o =  cdata.orange.crlbs(:,6).^0.5./cdata.orange.fits(:,6);
cdata.orange.ilv = cdata.orange.llv(:)./cdata.orange.fits(:,3);
cdata.orange.eps = abs(cdata.orange.fits(:,4)./cdata.orange.fits(:,5));
cdata.orange.snr =(cdata.orange.fits(:,3)./(cdata.orange.fits(:,3)+(2*6+1)^2*cdata.orange.fits(:,6)).^0.5);

% Apply Tolerance

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

fnames = fieldnames(cdata.red);
for i = 1:numel(fnames)
    cdata.red.(fnames{i})(ind,:) = [];
end
% Scream here for sanity

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

cdata.orange.xf = cdata.orange.xf(:);
cdata.orange.yf = cdata.orange.yf(:);
fnames = fieldnames(cdata.red);
for i = 1:numel(fnames)
    cdata.orange.(fnames{i})(ind,:) = [];
end
clear ind fnames
save(['toleranced\',fname(1:end-4),'_tol.mat']);
end