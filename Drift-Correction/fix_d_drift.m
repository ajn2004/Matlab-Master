% function func_spline_drift(fdname, pix_size, chunk_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Drift Correction on GPU
%
%  This is a script to address drift in FPALM images analyzed with Neural
%  Quhzx and other GPU oriented codes. The idea is that it will follow the
%  same math as the drift correction gui but
%
%
%
%
% AJN 7-21-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

pix_size = 50; % Final pixel size in nanometers
chunk_size = 1000; % number of frames to construct a partial render for correlation
mkdir('DC');

finfo = dir('*tol.mat');
if isempty(finfo)
    [fname, fpath] = uigetfile('*tol.mat');
    cd(fpath);
    finfo  = dir('*tol.mat');
end

for i = 1:numel(finfo)
<<<<<<< Updated upstream
=======
    try
>>>>>>> Stashed changes
    func_spline_drift(finfo(i).name,pix_size,chunk_size);
end