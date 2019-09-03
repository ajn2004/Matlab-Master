function zf = func_build_ramp(zs,framenumber,dz,range,j)
cyc = abs(range*1000/dz); % define cycle variable for looping
y = 0; %initialize y
cyc = abs(round(range*1000/dz));
x = 0:2*cyc;
zf = zs*0;
fnumber = framenumber + j;
y = 0;
dz = 0.8*dz;
clear dcoords
% close all
for i = 2:numel(x)
    y(i) = y(i-1)+dz*0.001;
    if mod(i-1,cyc) == 0
        dz = -dz;
    end
end

for i = 1:numel(zs(:,1))
    zf(i) = zs(i) + y(mod(fnumber(i),2*cyc)+1);
%     dcoords(i) = zf(i) + y(framenumber(i)+1)/q;
end