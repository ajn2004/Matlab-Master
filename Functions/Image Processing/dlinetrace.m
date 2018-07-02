function [path, totes] = dlinetrace()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Line Tracer
%
% A short script to return the line trace over several frames in hopes of
% being able to linearize and track motion along an axon
% AJN 6/28/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;


try
    c= scrub_config();
    i1 = readtiff()/em_gain(c.Gain);
catch lsterr
    i1 = readtiff();
end

[path, ~] = linetrace(sum(i1,3));

[~,~,o] = size(i1);
totes = [];
for i = 1:o
    [~,bright] = linetrace(i1(:,:,i),path);
    totes = [totes,bright(:)];
end
imagesc(totes)
colormap('jet')
title('Dynamic Line Trace Result')
xlabel('Framenumber')
ylabel('Path Index')