%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circular Density measure
%
% This script will measure the radial density function of regions returned
% by the gauss_density_2 program. It will return a group of radial density
% functions for all areas identified.
%
% THIS SCRIPT CONTAINS GUSTO!
%
% AJN 2/19/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all

% Factor of radius to look over
factor = 3;
param_thresh = 10;

% File Selection
[fname, fpath] = uigetfile('*regions.mat');
cd(fpath);
finfo = dir('*regions.mat');
figure('units','normalized','OuterPosition',[0 0 1 1])
index = [8, 19];
xq = 1:0.01:3;
for j = 1:numel(index)
    close all
    load(finfo(index(j)).name);
    % load(finfo(j).name);
    clear rtot dtot
    try
        cents = more_cents;
        radii = more_rads;
        areas = cat(1, s.Area);
        perims = cat(1, s.Perimeter);
        
        params = 4*pi*areas./perims.^2;
        
        rtot = [];
        dtot=[];
        
        for i = 1:numel(radii)
            
            if radii(i) > param_thresh
                [rs, ds] = circular_density(cents(i,1),cents(i,2),radii(i),im1,factor);
                
                rtot = [rtot;rs(:)./radii(i)];
                dtot = [dtot;ds(:)];
            end
            
        end
        
        figure('units','normalized','outerposition',[ 0 0 1 1])
        plot(rtot,dtot/(max(dtot(dtot<1000000))))
        title(finfo(index(j)).name)
        xlabel('Normalized Radius');
        ylabel('Normalized Density');
        % title(finfo(j).name)
    catch lsterr
    end
end
