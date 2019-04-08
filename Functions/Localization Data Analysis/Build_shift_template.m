% Sawtooth Grabber
% This script will grab a sawtooth wave from latex beed localization data
% to build a template for stage oscillating behavior.
clearvars;
close all;
clc
load('C:\Users\AJN Lab\Dropbox\Data\4-3-19 calibration\Analysis\latex_bead_dz10_r1_dast.mat');
cyc = 200;
offset = 3;
% zf = getdz(fits(:,4),fits(:,5),cal.z_cal);
zf = ncoords(:,3)*q;
tf = framenumber;

%segment data to get the average 'cycle'
x = tf(offset:offset+cyc);
zavg = zf(offset:offset+cyc);
count = 1;
while true
    try
        zavg = zavg + zf(count*cyc+offset:offset+cyc*(count+1));
        count = count +1;
    catch lsterr
        break
    end
end
zavg = zavg/count;
% plot(x,zavg,'.')
% hold on
plot(framenumber,zf,'.')
% hold on
func_build_shift_correct

a = polyfit(x(1:102),zavg(1:102),1)
a = polyfit(x(102:end),zavg(102:end),1)
save('dz10_r1_template.mat','offset','y')
