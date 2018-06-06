function ik = kfact(i1,varargin)
% this function is to implement the non-linear decomposition method called
% k-factor as described in Ilovitsh et al 2014 Biomedical Optics Express
% while k needs to be higher than a minimum value to ensure convergence,
% experiments have shown that choosing a value between 0.5 and 0.9 should
% reliably reach convergence of the geometric sum for the minimal values in
% the image

numvar = length(varargin);
if numvar > 3
    error('Too many variables');
elseif numvar < 3
    disp('Using Default Values of k = 0.9, M = 48 and h = 8');
    k = 0.9;
    M = 48;
    h = 8;
elseif numvar == 3
    k = varargin{1};
    M = varargin{2};
    h = varargin{3};
end

[m,n,o] = size(i1); % get image dimensions
ik = i1.*0; % preallocate output

% loop over all frames and apply the k-factor algorithm to each frame
for j = 1:o
    %linearize array
    isub = i1(:,:,j);
    fy = (isub(:)/max(isub(:))).'; % all values must be from 0 to 1
    
    for i = 1:h % k factor algorithm
        
        g(i,:) = fy >= (1+k^i)^-1; % determine the g value
        f(i,:) = (1 + k^i.*g(i,:))./(1+k^i); % determine the f value
        fy = fy./f(i,:); % reduce the remaining signal
    end
    rf = prod(f,1).';
    ik(:,:,j) = reshape(rf,m,n)*max(isub(:));
    clear rf fy isub g f
end