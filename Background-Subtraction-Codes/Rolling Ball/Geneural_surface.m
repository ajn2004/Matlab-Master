% Geneural Surface
%  A script for using a growing neural gas to learn the topology of data
%  sets. This is being specfically developed for use in 3D localization
%  data
close all; clearvars; clc;

%% User Defined Variables

% File Variables
% These variables control information regarding files to analyze

an_name = 'Geneural'; % leave empty to save file in data folder
f_start = 1;
f_end = 0;
graph = 'Y'; % set to 'Y' or 'y' if you want to see a final graph

% Data Variables
% these variables allow the user to select the data that they will be
% working with
% names will contain the strings of all variables to be loaded from the
% data set must be used as a names = {  ' name 1', 'name 2', ... ,' name
% n'} for n-dimensional data
names = {'ncoords','q','framenumber'};
dim = 3; % dimensionality of data

%% Neural Gas Variables
epsilon_b = 0.2;    % Scale factor to reduce 'step' of nearest node
epsilon_n = 0.006;  % scale factor to reduce 'step' of neighbor nodes
     amax = 50;     % Max age of a connection before cutting
     lamb = 100;    % cycle variable that chooses when to add nodes
    alpha = 0.5;    % Error reducer of newly created node
     decr = 0.995;  % Iteration error reduction
     params = [epsilon_b, epsilon_n, amax, lamb, alpha, decr];
     
     
% user Independent
[fname, fpath] = uigetfile('*.mat');
% addpath(pwd);
cd(fpath)
% % fname = 'well 1 cell 6_frames_1_10000_t_neuro_tol.mat';
% an_dir = [fpath,an_name,'\'];
mkdir(an_name);

finfo = dir('*.mat');
if f_end == 0
    f_end = numel(finfo);
end
an_dir = [fpath,an_dir];
for j = f_start:f_end
    close all
    for i = 1:numel(names)
        load(finfo(j).name,names{i});
    end
    r = str2num(fname(strfind(finfo(j).name,'_r')+2)); % get range of shift
    zf = func_shift_correct(ncoords(:,3)*q,framenumber,r).'; % shift correct
    data = [ncoords(:,1)*q,ncoords(:,2)*q,zf(:)]; % create 'data' variable
    [w, T, A] = func_gng(data,params)