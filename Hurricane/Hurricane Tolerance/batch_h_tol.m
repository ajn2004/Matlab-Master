%Batch H_tol
clearvars; close all; clc;
mkdir('toleranced')
files = dir('*dast*');

for i = 1:numel(files)
    func_batch_h_tol(files(i).name);
end
clear files i