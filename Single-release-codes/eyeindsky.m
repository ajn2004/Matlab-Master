function [cents, thrsh] = eyeindsky(i1,molish)
% this is a function to find molecules in noisey images by looking at the
% ratio of a 'signal' image and the 'background' image

cents = []; % preallocate centers variable
rat1 = bandpass(i1(:,:,2)./i1(:,:,1)); % get the ratio and apply a spatial bandpass filter

% kill border
rat1(1:6,:) = 0;
rat1(:,1:6) = 0;
rat1(end-5:end,:) = 0;
rat1(:,end-5:end) = 0;

sr1 = sort(rat1(:)); % sort pixels by ascending order)
numpix = 2*molish*9*9; % molish is the estimated number of molecules, which are hypothesized to give the brightest pixels
                       % here we will multiply that by 2 for good measure
                       
thresh = sr1(end-numpix+1); % the threshold value will be the value of the pixel numpix away from maximum

dps = get_das_peaks(rat1,thresh); % get the local maxima in the ratio image

rps = dps.*i1(:,:,2);  % multiplying the peaks image by the ratio image gives a sparse image with only local maxima values
sr2 = sort(rps(:));  % repeat sort
% imagesc(rat1);
thrsh = sr2(end-molish+1); % select only the molish number of regions

[row, col] = find(rps >= thrsh);
cents = [col, row];