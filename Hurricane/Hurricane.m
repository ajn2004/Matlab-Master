close all; clearvars; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hurricane
%
% Hurricane will be a full localization software suite by the time I am
% done. Right now it is a batch program that allows the user to run
% Da_Storm over every tiff in a folder. This is done similar to the way the
% "Bracket" localization algorithms perform
%
% AJN 7/25/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Hurricane(thresh)

%% User Controlled Area
emg = 0; % Enter the EM Gain value used on the camera
pix2pho = em_gain(emg);    %Pixel to photon ratio
q = 0.122;          % Pixel size in um
pixw = 6;       % radius to localize (final image size is (2*pixw+1)^2 pixels)
an_dir = 'Analysis'; % name of analysis directory
angle = 0; %approximate astig rotation in degrees
sv_im = 'n'; % set to y/Y to save image of localizations
thresh = 15;

%% Optionals
  % This section is dedicated to a list of variables for the user to select
  % 1 indicates go 0 indicates do not
  savewaves = 0;
  showlocs = 0;
  savepsfs = 0;
  saverb = 0;
  two_color = 1;
  varys = [savewaves, showlocs, savepsfs, saverb, two_color];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% END USER CONTROL JUST RUN IT AND SELECT A FILE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
p = mfilename('fullpath');
[fpath, fname, fext] = fileparts([p, '.m']);
addpath(fpath);
addpath([fpath,'\da_c']);
addpath([fpath,'\da_functions']);

% [dname, dpath] = uigetfile('*.tif');
% cd(dpath); % we change to the path we are going to use
dpath = pwd;

try
load('back_subtract.mat');
catch lsterr
    mi1 = 0;
end

files = dir('*.tif');
% if ~isempty([fpath,'\',an_dir])
% sendit2(an_dir);
% end
mkdir(an_dir);
%% Localize the files with the thresholds found
% thresh = findathresh(files,pix2pho,mi1);
% mi1 = 0;
if varys(1) == 1
    mkdir('Waves');
elseif varys(3) == 1
    mkdir('psfs');
elseif varys(4) == 1
    mkdir('Rolling_Ball');
end
for i = 1:numel(files)
    tic
    func_da_storm(files(i).name, dpath, an_dir, q, pix2pho, pixw,thresh, angle, sv_im, mi1, varys);
    clc;
    disp(['File number ' , num2str(i) , ' out of ', num2str(numel(files))]);
    t(i) = toc;
    ajn_wait(t,i,numel(files));
    close all
end



% end

