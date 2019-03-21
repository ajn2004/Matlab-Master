function func_batch_h_tol(fname)
% batch routine for files
load('Tolfile.mat');
load(fname);

%Define fractionals
fr_N =  crlbs(:,3).^0.5./fits(:,3);
fr_sx = crlbs(:,4).^0.5./fits(:,4);
fr_sy = crlbs(:,5).^0.5./fits(:,5);
fr_o =  crlbs(:,6).^0.5./fits(:,6);
ilv = llv(:)./fits(:,3);
eps = abs(fits(:,4)./fits(:,5));

% Apply Tolerances
ind = fits(:,3) > N_tol(1) & fits(:,3) < N_tol(2); % Photon Tolerance
ind = ind & abs(ncoords(:,3)) < absz;
ind = ind & fits(:,4) > s_tol(1) & fits(:,4) < s_tol(2); % Photon Tolerance
ind = ind &fits(:,5) > s_tol(1) & fits(:,5) < s_tol(2); % Photon Tolerance
ind = ind & q*crlbs(:,1).^.5 < lat_max & q*crlbs(:,2).^.5 < lat_max; % Lateral Uncertainty Tolerance
ind = ind & ilv > iln; % llv tolerance
ind = ind & fr_N < frac_lim & abs(fr_o) < off_frac; % Fraction photon tolerance
ind = ind & fr_sx < frac_lim & fr_sy < frac_lim; % Fraction width tolerance
ind = ind & eps <= maxe & eps >= mine;

fits(logical(1-ind),:) = [];
crlbs(logical(1-ind),:) = [];
llv(logical(1-ind)) = [];
ncoords(logical(1-ind),:) = [];
framenumber(logical(1-ind)) = [];
clear ind N_tol s_tol iln frac_lim lat_max fr_N fr_sx fr_sy fr_o ilv absz
save(['toleranced\',fname(1:end-4),'_tol.mat']);
end