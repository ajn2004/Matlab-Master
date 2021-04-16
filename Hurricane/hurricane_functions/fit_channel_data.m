function cdata = fit_channel_data(varargin)

% populate variables from varargin
iloc = varargin{1}; %Images to localize
fnum = varargin{2}; % Array of framenumber
cents = varargin{3}; % Array of central pixels
cal = varargin{4}; % Calibration data structure
color  = varargin{5}; % string variable for data structure
q = cal.q;
if nargin == 6 % If 6 arguments then we have a prepopulated cdata variable
    cdata = varargin{6};
end
% A funciton to analyze channel data based off the 'color' variable value
% iloc is a list of localization images, fnum is the array of frames each
% loc was found, cents is an array of central positions, cal is the
% calibration data structure, color is a string variable 'red' or 'orange'
[fits, crlbs, llv, framenumber, psfs] = slim_locs(iloc, fnum, cents, cal.(color).ang);
            
% As everywhere in the equations used sigma is squared, we can without
% loss of generality make these fits positive definite
fits(:,4) = abs(fits(:,4));
fits(:,5) = abs(fits(:,5));
[m,n,o] = size(iloc);
% Z calculations
zf = get_spline_z(fits(:,4),fits(:,5),cal.(color)); % New z_registration based off spline 3d calibration
% remove failed z assignments
ind = zf == -5;
fits(ind,:) = [];
crlbs(ind,:) = [];
llv(ind) = [];
framenumber(ind) = [];
zf(ind) = [];
psfs(:,ind) = [];
o = numel(psfs(1,:));
psfs = reshape(psfs,m,n,o);

% Put data into cdata structure
for i = 1:6
    cdata.(color).fits(:,i) = fits(:,i);
    cdata.(color).crlbs(:,i) = crlbs(:,i);
end
cdata.(color).llv = llv;
cdata.(color).framenumber = framenumber;
cdata.(color).psfs = psfs;

ncoords = make_astigmatism_corrections([cdata.(color).fits(:,1:2),zf/q],cal.(color),cal.q);
% Assign fixed coordinates which are in microns at this point
cdata.(color).xf = ncoords(:,1);
cdata.(color).yf = ncoords(:,2);
cdata.(color).zf = ncoords(:,3);

index = cdata.(color).xf  < 0 | cdata.(color).yf <0 ;
field_names = fieldnames(cdata.(color));

if sum(index)>0
    cdata.(color).fits(index,:) = [];
    cdata.(color).crlbs(index,:) = [];
end

for k=3:numel(field_names)
    if sum(index) > 0
        cdata.(color).(field_names{k})(index) = [];
    end
    
end
sx_index = cal.(color).z0s > -0.5 & cal.(color).z0s < 0.5;
sx_data = [cal.(color).sx(sx_index);cal.(color).sy(sx_index)].';
nn_index = knnsearch(sx_data,[cdata.(color).fits(:,4),cdata.(color).fits(:,5)]);
for i = 1:numel(cdata.(color).xf)
    cdata.(color).distance(i,1) =  ((cdata.(color).fits(i,4) - sx_data(nn_index(i),1)).^2 + (cdata.(color).fits(i,5) - sx_data(nn_index(i),2)).^2)^0.5;
end
end