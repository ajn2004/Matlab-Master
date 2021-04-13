function cdata_out = apply_tolerance_to_cdata(cdata_in, tol, color)
cdata_out = cdata_in;
%Eliminate out of view molecue
fnames = fieldnames(cdata_in.(color));
ind = cdata_in.(color).xf < 0 | cdata_in.(color).xf > 200 |cdata_in.(color).yf < 0 | cdata_in.(color).yf > 200;
for i = 1:numel(fnames)
    cdata_in.(color).(fnames{i})(ind,:) = [];
end
%Define fractionals
tol.r.dist = 0.5;
cdata_in.(color).fr_N =  cdata_in.(color).crlbs(:,3).^0.5./cdata_in.(color).fits(:,3);
cdata_in.(color).fr_sx = cdata_in.(color).crlbs(:,4).^0.5./cdata_in.(color).fits(:,4);
cdata_in.(color).fr_sy = cdata_in.(color).crlbs(:,5).^0.5./cdata_in.(color).fits(:,5);
cdata_in.(color).fr_o =  cdata_in.(color).crlbs(:,6).^0.5./cdata_in.(color).fits(:,6);
cdata_in.(color).ilv = cdata_in.(color).llv(:)./cdata_in.(color).fits(:,3);
cdata_in.(color).eps = abs(cdata_in.(color).fits(:,4)./cdata_in.(color).fits(:,5));
cdata_in.(color).snr =(cdata_in.(color).fits(:,3)./(cdata_in.(color).fits(:,3)+(2*6+1)^2*cdata_in.(color).fits(:,6)).^0.5);
%% Apply  red Tolerance
 cdata_in.(color).framenumber = cdata_in.(color).framenumber(:);
% Build Red
ind = cdata_in.(color).fits(:,3) > tol.r.N_tol(1) & cdata_in.(color).fits(:,3) < tol.r.N_tol(2); % Photon Tolerance
ind = ind & cdata_in.(color).snr >= tol.r.minsnr; % Signal to Noise tolerances
ind = ind & abs(cdata_in.(color).framenumber - mean(tol.r.flims)) <= diff(tol.r.flims)/2; % Framenumber tolerances
ind = ind & abs(cdata_in.(color).zf-mean(tol.r.zlims)) <= diff(tol.r.zlims)/2; % axial tolerances
ind = ind & abs(cdata_in.(color).fits(:,6)-mean(tol.r.offlim)) <= diff(tol.r.offlim)/2;   % offset tolerances
ind = ind & (abs(cdata_in.(color).fits(:,4)).*abs(cdata_in.(color).fits(:,5))).^0.5 > tol.r.s_tol(1) & (abs(cdata_in.(color).fits(:,4)).*abs(cdata_in.(color).  fits(:,5))).^0.5 < tol.r.s_tol(2); % sigma Tolerance
ind = ind & tol.q*cdata_in.(color).crlbs(:,1).^.5 < tol.r.lat_max & tol.q*cdata_in.(color).crlbs(:,2).^.5 < tol.r.lat_max; % Lateral Uncertainty Tolerance
ind = ind & cdata_in.(color).ilv > tol.r.iln; % llv tolerance
ind = ind & cdata_in.(color).fr_N < tol.r.frac_lim & abs(cdata_in.(color).fr_o) < tol.r.off_frac; % Fraction photon tolerance
ind = ind & cdata_in.(color).fr_sx < tol.r.frac_lim & cdata_in.(color).fr_sy < tol.r.frac_lim; % Fraction width tolerance
ind = ind & cdata_in.(color).distance < tol.r.dist;
rind = ind; % Index represents all molecules that PASS tolerances
fnames = fieldnames(cdata_in.(color));
ind = logical(1-rind); % to remove failed localizations we want 1-pass
for i = 1:numel(fnames) % remove data points
    cdata_in.(color).(fnames{i})(ind,:) = [];
end
 cdata_out = cdata_in;
