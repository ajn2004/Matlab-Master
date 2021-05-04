% The Perfect Storm
% This is a batch script that will fully analyze a collection of folders
% specified by the user. We will be calling the batch version of all our
% functions so the user variables setting will be a bit expansive
% AJN 9-27-2020
clearvars;
close all;
clc;

%% Regular Change User Variables
folder_names_to_analyze = {'5-3-21 gpi-vglutmeos'};
align_color = 'red';

%% Set and forget Variables
% Hurricane Variables
emg = 0; % Enter the EM Gain value used on the camera
pix2pho = em_gain(emg);    %Pixel to photon ratio
q = 0.122;          % Pixel size in um
pixw = 6;       % radius to localize (final image size is (2*pixw+1)^2 pixels)
an_dir = 'Analysis'; % name of analysis directory
angle = 0; %approximate astig rotation in degrees
sv_im = 'n'; % set to y/Y to save image of localizations
framechunk_drift_correct = 1000;
thresh = 10;

%% Perform the Hurricane Process
% Hurricane Optionals
% This section is dedicated to a list of variables for the user to select
% 1 indicates go 0 indicates do not
savewaves = 0;
showlocs = 0;
savepsfs = 0;
saverb = 0;
two_color = 1; % two color code is as follows 1 = 2 color (algorithm decides order), 2= red only 3 = orange only 0 = no frame blocking
varys = [savewaves, showlocs, savepsfs, saverb, two_color];

for l = 1:numel(folder_names_to_analyze)
try
    folder_to_analyze = ['G:\Dropbox\Data\', folder_names_to_analyze{l}, '\'];
    cd(folder_to_analyze);
    
catch lsterr
    folder_to_analyze = ['D:\Dropbox\Data\', folder_names_to_analyze{l}, '\'];
    cd(folder_to_analyze);
end
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
% Normalize Name information
files = rename_problem_files(files);
if varys(1) == 1
    mkdir('Waves');
elseif varys(3) == 1
    mkdir('psfs');
elseif varys(4) == 1
    mkdir('Rolling_Ball');
end
disp('Hurricane')

lost_inds = [];
mkdir('Analysis\Raw\');
mkdir('Analysis\Tol\');
mkdir('Analysis\Fin\');
for i = 1:numel(files)
    if isempty(strfind(files(i).name,'scan'))
        try
            filename = [folder_to_analyze,files(i).name];
            disp(filename)
            func_da_storm_ps(files(i).name, q, pix2pho, pixw,thresh, angle, sv_im, mi1, varys);
        catch lsterr
            disp(lsterr.message)
            lost_inds = [lost_inds;i];
        end
    end
end
disp('Tolerance')
% 
% % Batch H_tol
% cd(folders_to_analyze{l});
toleranced_files = dir('Analysis\raw\*dast.mat');

disp('files')
disp(lost_inds)
disp('Were lost during hurricane')
lost_inds = [];
% Batch Scan Correction
folder = [folder_to_analyze,'Analysis\Raw\'];
for i = 1:numel(toleranced_files)
    if isempty(strfind(toleranced_files(i).name,'scan'))
        try
            
            file = toleranced_files(i).name;
%         filename = [folders_to_analyze{l},'Analysis\',files(i).name];
        func_batch_h2_tol_ps(folder, file);
%         delete(filename)
        catch lsterr
            disp(lsterr)
            lost_inds = [lost_inds;i];
        end
    end
end

disp('files')
disp(lost_inds)
disp('Were lost during tolerance')
cd('Analysis\');
folder = [folder_to_analyze,'Analysis\Tol\'];
drift_and_cluster_correct_tol_folder(folder,'red',1000);

% image_path = folders_to_analyze{l};
% lost_inds = [];
% for i = 1:numel(files)
%     if isempty(strfind(files(i).name,'scan'))
%         try
%             image_file_name = [image_path, files(i).name];
%             image_ruler_name = [image_path, files(i).name(1:end-8),'scan.tif'];
%             filename = [folders_to_analyze{l},'Analysis\',files(i).name(1:end-9),'_dast_tol.mat'];
%             file_list = {filename,image_file_name,image_ruler_name};
%             t(i) = laser_scan_correction_avg_ps(file_list);
% %             t(i) = laser_scan_correction_ps(file_list);
% %             delete(filename)
%         catch
%             lost_inds = [lost_inds;i];
%         end
%     end
% end
% disp('files')
% disp(lost_inds)
% disp('Were lost during scan correction')
% lost_inds = [];
% for i  = 1:numel(files)
%     if isempty(strfind(files(i).name,'scan'))
%         filename = [folders_to_analyze{l},'Analysis\',files(i).name(1:end-9),'_dast_tol_sc.mat'];
%         try
%             load(filename);
%             
%             try
%                 cdata.orange = converge_duplicates(cdata.orange,2, 1);
%             catch
%             end
%             try
%                 cdata.red = converge_duplicates(cdata.red,2, 0);
%             catch
%             end
%             save([folders_to_analyze{l},'Analysis\',files(i).name(1:end-4),'_full.mat'],'cdata','tol','cal');
%             delete(filename)
%         catch
%             
%         end
%     end
% end
% disp('files')
% disp(lost_inds)
% disp('Were lost during convergance')
end
