function check_tag_load(folder,file_name)
% This will become the fundamental way to load localization data
% load file
load([folder,'\',file_name]);
index_file_path = '..\..\Composites\';
index_file_name = 'data_index.mat';
load([index_file_path,index_file_name])
% Check if file contains an id listing, if it does not we will add it to
% the repository
if ~exist('file_id_number','var') % Loop executes in absence of a number
    
    % New ID number is 1 + previously assigned
    file_id_number = numel(data_index) + 1;
    % Acknowledge new file to user
    disp(['New file detected'])
    % Save file id w/ file
    save([folder,'\',file_name],'cdata','cal','tol','file_id_number')
    % Record file name in proper location in data index object
    data_index(file_id_number).name = file_name;
    % Record file location in data index object
    data_index(file_id_number).folder = folder;
end
% Print File information
disp(['Loading file ', file_name, ' from folder ', folder]);
disp(['ID is ', num2str(file_id_number)]);
% save updated data index object
save([index_file_path,index_file_name],'data_index');
% Load data variables into proper workplace scope
get_coords_in_workspace(cdata)
assignin('base','cal',cal);
assignin('base','cdata',cdata);
assignin('base','data_index',data_index);
assignin('base','file_id_number',file_id_number);
end