function check_tag(folflder,file_name)
% This will become the fundamental way to load localization data
% load file
load([folflder,'\',file_name]);
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
    save([folflder,'\',file_name],'cdata','cal','tol','file_id_number')
    % Record file name in proper location in data index object
    data_index(file_id_number).name = file_name;
    % Record file location in data index object
    data_index(file_id_number).folder = folflder;
end

% save updated data index object
save([index_file_path,index_file_name],'data_index');

end