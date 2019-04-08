function [m,n] = hd_stats(res)
% A quick look up table for HD pixel sizes
% AJN 4-8-19

switch res
    case 720
        m = 1280;
        n = 720;
    case 1080
        m = 1920;
        n = 1080;
    case 4000
        m = 3840;
        n = 2160;
    case 4001
        m = 4096;
        n = 2160;
end