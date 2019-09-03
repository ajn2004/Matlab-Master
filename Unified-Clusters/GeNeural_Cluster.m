%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GeNeural Cluster
%
%   This script utilizing functions such as neural gas and
%   to accomplish object identification in a user specified file
%
%
%   AJN 5/11/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% User defined Variables

% File Variables
% These variables control information regarding files to analyze
mkdir('Geneural')
an_dir = 'Geneural\'; % leave empty to save file in data folder
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

% Neural Gas Variables
% Resonable values for each variable are provided in comments

nodes = 500; % The number of nodes making up your gas  reasonable: 200
epsilon_i = 0.3; % step size in weight vector learning eq 1. should be between 0 and 1
epsilon_f = 0.01;  % reasonable values are e_i = 0.3 and e_f = 0.01
% This variable controls the extent of the step size,
% too large inital variable will result into
% unreasonable clumping of many of the nodes

lambdai = 250;  % decay rate vector for exponential decay distance updating
lambdaf = 0.5;  % This variable affects how much the nodes will expand into the less dense regions
% If the initial value is too low relative to your number
% of nodes, you will see more nodes left outside the region
% of interest and not connected at all
% Too high final values will not allow for proper filling of
% clustered strutcure. Resonable values lambdai = 0.5*nodes,
% lambdaf = 0.5

lifeti = 20;  % maximum lifetime of any connection without being refreshed
lifetf = 80;   % This variable controls how connected the web is in the end, if these values are too high
% You will notice each node having too many connections, if
% too low you will not have connected structure
% reasonable values are lifeti = 20 lifetf = 80

tmax = 4000;   % this is the final iteration  reasonable: 40000
ng_params0 = [epsilon_i,epsilon_f, lambdai,lambdaf, lifeti, lifetf, nodes];
% Analysis Variables

dr = 10;
rmax = 500;
% user Independent
[fname, fpath] = uigetfile('*.mat');
% addpath(pwd);
% cd(fpath)
% % fname = 'well 1 cell 6_frames_1_10000_t_neuro_tol.mat';

finfo = dir('*.mat');
if f_end == 0
    f_end = numel(finfo);
end
for j = f_start:f_end
    close all
    for i = 1:numel(names)
        load(finfo(j).name,names{i});
    end
    %% User Specified Data preparation
    % this section allows the user to specify the process of creating a
    % data matrix which must take the form of data(datapoint,
    % datadimension)
    r = str2num(fname(strfind(finfo(j).name,'_r')+2));
    zf = func_shift_correct(ncoords(:,3)*q,framenumber,r).';
    data = [ncoords(:,1)*q, ncoords(:,2)*q,zf];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Thanks user, you're shouldn't need to change anything else
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% Neural Gas Section
    set(0,'RecursionLimit',5000)
    % Run the neural gas algorithm
%     [w, C, T] = func_neural_gas(data, nodes, epsilon_i, epsilon_f, lambdai, lambdaf, lifeti, lifetf, tmax, graph);
    % w is a matrix containing the 'positions' of the neural nodes
    % C is a matrix showing connectivity between node i and node j at
    % C(i,j) where C(i,j) is 1 if i and j are connected and 0 if not
    % T is an age map that shows the age of the connection between nodes i
    % and j at T(i,j) if C(i,j) == 1
    
    %% Structural Analysis
    
    % Radial Distribution analysis
    [R, g] = radial_df(w*1000, dr, rmax);
    figure
    subplot(2,2,1);plot(R,mean(g,2))
    title('Pair Distribution Curve for Selecting Clustering');
    xlabel('Radial distance in nm');
    ylabel('g(r)');
    
    %     % Web Cutting
    %     [Cnew, wnew, clus_id] = web_cutter(w,C, data, [1, dr, rmax], R);
    
    % %     %% Data Saving
    M = g;
    varys = {'w','C','T'}; %variables to be saved in addtion to data variables
    savename1 = [an_dir,finfo(j).name(1:end-4),'_NG'];
    eval(quicksave(savename1,'.mat',names,varys));
    %
end