%% JImage splitter
% This script is to separate data in Lorena's files from widefield to
% single molecule
clearvars; close all; clc;


[fpath, fname, fext] = fileparts([pwd, '.m']);
files = dir('*tif');
mkdir('split');

for i = 1:numel(files)
    i1 = readtiff(files(i).name);
    [m,n,o] = size(i1);
    index = 2:2:o;
    writetiff(i1(:,:,index),[fpath,'\Lorena and Andrew\split\',files(i).name(1:end-4),'_evens']);
    writetiff(i1(:,:,index-1),[fpath,'\Lorena and Andrew\split\',files(i).name(1:end-4),'_odds']);
end