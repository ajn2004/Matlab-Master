close all; clear all; clc;
noc = 3;
load('I:\Data\9-10-15 Mito TCS FPALM\Analysis\Toleranced\Combined\Drift Corrected\Mito Analysis\Control_results_wout_3_and_8.mat');
areatot = Area;
perimtot = Perimeter;
eccen = Eccentricity;
ediam = EquivDiameter;
farea = FilledArea;
soli = Solidity;
m = numel(Area);
load('I:\Data\9-10-15 Mito TCS FPALM\Analysis\Toleranced\Combined\Drift Corrected\Mito Analysis\Treatment_results.mat');

areatot = vertcat(areatot, Area);
perimtot = vertcat(perimtot,Perimeter);
eccen = vertcat(eccen,Eccentricity);
ediam = vertcat(ediam,EquivDiameter);
farea = vertcat(farea,FilledArea);
soli = vertcat(soli,Solidity);

data = [areatot./max(areatot), perimtot./max(perimtot), eccen, ediam./max(ediam), farea./max(farea), soli];

[clus_coords, index, M] = funkmeanscluster(noc, data, 'y');


plot(index, '.b');
hold on
plot([m,m],[1,3], 'r')
title('Cluster identification by index');
xlabel('Index of mitochondria');
ylabel('Cluster assignment');
hold off