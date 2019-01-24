% Data Combiner
% Combines large data sets into a single data set to be used on toleranced
% files
%
% 9/17/15 AJN
try
while true
clear all
close all
clc

frames = 4085;
nfiles = 1;
mkdir('tot');

% Use a try catch loop to load as many files as necessary, handles variable
% number of files
try
    while true
        
        [fname_temp, fpath_temp] = uigetfile('*dast*','Select a file in chronological order, cancel to continue');
        cd(fpath_temp);
        finame{nfiles} = cellstr(fname_temp);
        fpath{nfiles} = cellstr(fpath_temp);
        nfiles = nfiles +1;
    end
catch lasterr
    nfiles=nfiles - 1;
    disp(['Combining ', num2str(nfiles), ' files']);
end

% preallocate initial arrays
tcoords =[];
tcrlbs = [];
tfits = [];
tfnumb = [];
tlv = [];

imfiles = [];
total_mol = 0;
%loop over all files combining relevant localization data
for jk = 1:nfiles
    load([char(fpath{jk}), char(finame{jk})]);
    if jk == 1
        tot_name = char(finame{jk}); 
    end
    frames_last_file(jk+1) = framenumber(end);
    tcoords = [tcoords;ncoords];
    tcrlbs = [tcrlbs;crlbs];
    tfits = [tfits;fits];
    tfnumb = [tfnumb; framenumber+sum(frames_last_file(1:jk))];
    tlv = [tlv; llv];
    imfiles{jk} = char(finame{jk});
    clear fits crlbs fname framenumber llv ncoords
end
total_mol = numel(tlv);
ncoords = tcoords;
crlbs = tcrlbs;
fits = tfits;
framenumber = tfnumb;
llv = tlv;
clear tcoords tcrlbs tfits tfnumb tlv

save([char(fpath{jk}),'\tot\',tot_name(1:end-4),'_tot.mat'])
end
catch lsterr
end