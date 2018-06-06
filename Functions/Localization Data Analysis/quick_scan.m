% quick_scan
% Quickly scan through localization files to find the one you like

while true
    [fname, fpath] = uigetfile('*.mat');
    load(fname);
    plot(xf_all,yf_all,'.');
    axis equal
    title(fname);
    waitforbuttonpress;
end