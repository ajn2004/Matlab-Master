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
an_dir = ['F:\Data\2-17-16 live cell mito tcs\Analysis\Tol\Clusters\']; % leave empty to save file in data folder
f_start = 1;
f_end = 1;
graph = 'Y'; % set to 'Y' or 'y' if you want to see a final graph


% Data Variables
% these variables allow the user to select the data that they will be
% working with
% names will contain the strings of all variables to be loaded from the
% data set must be used as a names = {  ' name 1', 'name 2', ... ,' name
% n'} for n-dimensional data
names = {'xf_all', 'yf_all','q'};
dim = 2; % dimensionality of data

% Neural Gas Variables
% Resonable values for each variable are provided in comments

nodes = 1600; % The number of nodes making up your gas  reasonable: 800
epsilon_i = 0.3; % step size in weight vector learning eq 1. should be between 0 and 1
epsilon_f = 0.08;  % reasonable values are e_i = 0.3 and e_f = 0.01
                    % This variable controls the extent of the step size,
                    % too large inital variable will result into
                    % unreasonable clumping of many of the nodes

lambdai = 500;  % decay rate vector for exponential decay distance updating
lambdaf = 0.5;  % This variable affects how much the nodes will expand into the less dense regions
                % If the initial value is too low relative to your number
                % of nodes, you will see more nodes left outside the region
                % of interest and not connected at all
                % Too high final values will not allow for proper filling of
                % clustered strutcure. Resonable values lambdai = 80,
                % lambdaf = 0.5

lifeti = 20;  % maximum lifetime of any connection without being refreshed
lifetf = 40;   % This variable controls how connected the web is in the end, if these values are too high
             % You will notice each node having too many connections, if
             % too low you will not have connected structure
               % reasonable values are lifeti = 20 lifetf = 80

tmax = 50000;   % this is the final iteration  reasonable: 40000

% Analysis Variables
dr = 20;
rmax = 4000;
% user Independent
[fname, fpath] = uigetfile('*tol.mat');
addpath(pwd);
cd(fpath)

finfo = dir('*tol.mat');
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
    data = [xf_all*q,yf_all*q];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Thanks user, you're shouldn't need to change anything else
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Neural Gas Section
    % Run the neural gas algorithm
    [w, C, T] = func_neural_gas(data, nodes, epsilon_i, epsilon_f, lambdai, lambdaf, lifeti, lifetf, tmax, graph);
    % w is a matrix containing the 'positions' of the neural nodes
    % C is a matrix showing connectivity between node i and node j at
    % C(i,j) where C(i,j) is 1 if i and j are connected and 0 if not
    % T is an age map that shows the age of the connection between nodes i
    % and j at T(i,j) if C(i,j) == 1
    
    %% Structural Analysis
    
    % Radial Distribution analysis
    [R, g] = radial_df(w*1000, dr, rmax);
    figure
    plot(R,mean(g,2))
    title('Pair Distribution Curve for Selecting Clustering');
    xlabel('Radial distance in nm');
    ylabel('g(r)');
    
    
    
    % Nearest Neighbor Analysis;
    % Runs analysis returning the nearest neighbor distance for the ith
    % node
    
    nn_dist = func_nn_nodes(nodes,w);
    [counts, edges] = histcounts(nn_dist);
    mids = diff(edges)/2 + edges(1:end-1);
    
    
    % find NN_peak and widths for and threshold
    xind = find(counts == max(counts));
    beta0(1) = mids(xind);
    beta0(2) = counts(xind);
    sigind = find(counts(xind:end) <= max(counts) * 0.1353,1) + xind - 1;
    beta0(3) = (mids(sigind) - beta0(1))/2;
    beta0(4) = mean(counts(sigind:end));
    [betafit,resid,J,COVB,mse] = nlinfit(mids,counts.',@gauss1d,beta0);
    nn_max = 0.241;
    
    % remove under connected nodes
    indexout = zeros(nodes,1);
    for i = 1:nodes
        if sum(C(i,:)) < 2
            indexout(i) = 1;
        end
    end
    index = nn_dist <nn_max & logical(1-indexout);
    
    hold on
    plot(w(index,1),w(index,2),'b.', 'markersize',30)
    hold off
    drawnow
    
    %% Data Saving
    varys = {'w','C','T','nn_dist','index'}; %variables to be saved in addtion to data variables
    savename1 = [an_dir,finfo(j).name(1:end-4),'_clustered'];
%     eval(quicksave(savename1,'.mat',names,varys));
    
end