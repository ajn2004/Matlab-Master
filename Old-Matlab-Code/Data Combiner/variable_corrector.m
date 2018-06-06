% variable converter
close all
clear all
clc

fpath = 'I:\Data\10-15-15 Live Cell Mito FPALM\Analysis\Tol combined\';
cd(fpath);
files = dir('*.mat');

for kj = 1 : numel(files)
    load(files(kj).name);

        total_molecules = numel(yf_all);

     
    save(files(kj).name);
    clearvars -except fpath files kj
    kj
end
