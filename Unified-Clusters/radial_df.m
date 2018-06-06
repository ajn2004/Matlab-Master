function [R, g] = radial_df(X, dr, rmax)

    R = 1:dr:rmax; % looping over 3nm to 1um
    
   
    
    % convert spatial variables to nanometers
    xf_nm = X;
    yf_nm = X;
    cents_nm = X;
    
    %determine mean density
    minx = min(xf_nm);
    maxx = max(xf_nm);
    miny = min(yf_nm);
    maxy = max(yf_nm);
    
    maxr = mean([(maxx-minx)/2, (maxy-miny)/2]);
    rho = numel(xf_nm)/(pi*maxr^2);
    g = zeros(numel(R),numel(cents_nm(:,1)));
    % loop over all center positions
    figure
    for i = 1:numel(cents_nm(:,1))
        i
        rs = ((xf_nm - cents_nm(i,1)).^2 + (yf_nm - cents_nm(i,2)).^2).^0.5; %change coords to be radially away from point of interest
        count = 1;
        % radial Distribution calculation
        for j = 1:numel(R)
            index = find(rs > R(j) & rs< R(j)+ dr); % Find localization in shell of dR located R(j) away from point of interest
            n = numel(index);
            %         if j - dr/2 < 0
            %             sa = rho * pi * (j + dr/2)^2;
            %         else
            sa = rho*pi*((R(j)+ dr)^2-(R(j))^2); % normalizing factor
            %         end
            g(j,i) = n/sa; % RDF value at R(j)
            count = count+1;
        end
    end
