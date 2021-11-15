function smooth_surface = smooth_surface_data(varargin)

% Handle optional input variables
if nargin < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
end
% parse cdata out from optional inputs
[args,pvpairs] = parseparams(varargin);
surface_data = args{1};

% Check pvpairs for optional variables and assign them
if sum(strcmpi(pvpairs(1:end-1),'max_dist'))
    ind = find(strcmpi(pvpairs(1:end-1),'max_dist'));
    max_dist = pvpairs{ind+1};
end

if sum(strcmpi(pvpairs(1:end-1),'k'))
    ind = find(strcmpi(pvpairs(1:end-1),'k'));
    k = pvpairs{ind+1};
end

% build any missing variables
if ~exist('k','var')
    k = 4;
end
if ~exist('max_dist','var')
    max_dist = 0.07;
end

smooth_surface = [];
number = [];
[idx, d] = knnsearch(surface_data, surface_data,'k',k);
for i = 1:numel(surface_data(:,1))
    % Loop over all surface molecules
%     index = ((surface_data(:,1) - surface_data(i,1)).^2 + (surface_data(:,2) - surface_data(i,2)).^2 + (surface_data(:,3) - surface_data(i,3)).^2).^0.5 <= max_dist;
    smooth_surface(i,:) =  mean(surface_data(idx(i,:),:),1);
%     if sum(index) > k
%         smooth_surface =  [smooth_surface;mean(surface_data(index,:))];
%         number = [number;sum(index)];
%     end
end