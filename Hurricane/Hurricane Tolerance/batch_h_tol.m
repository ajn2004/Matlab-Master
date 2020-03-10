%Batch H_tol
clearvars; close all; clc;
mkdir('toleranced')
files = dir('*dast*');

for i = 1:numel(files)
    try
        load(files(i).name,'cdata');
        x = cdata;
        func_batch_h2_tol(files(i).name);
    catch lsterr
        func_batch_h_tol(files(i).name);
    end
end
clear files i