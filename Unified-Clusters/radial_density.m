%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radial distribution Function
%
% Given centers this will perform a radial distribution calculation
% AJN  2/26/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

[fname, fpath] = uigetfile('*clustered*.mat');
cd(fpath);
finfo = dir('*clustered*.mat');

for kj = 1:numel(finfo)
    load(finfo(kj).name);
    % Choose value over which R will scan
    dr = 10; % the shell of the dR
    R = 1:dr:4000; % looping over 3nm to 1um
    
    
    xf_nm = xf_all*q * 1000;
    yf_nm = yf_all*q * 1000;
    cents_nm = w(index,:)*1000;
    
    %determine mean density
    minx = min(xf_nm);
    maxx = max(xf_nm);
    miny = min(yf_nm);
    maxy = max(yf_nm);
    
    maxr = mean([(maxx-minx)/2, (maxy-miny)/2]);
    rho = numel(xf_nm)/(pi*maxr^2);
    g = zeros(numel(R),numel(cents_nm(:,1)));
    % loop over all center positions
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
    savename = [finfo(kj).name(1:end-10),'_rdf.mat'];
    save(savename,'g','R')
end