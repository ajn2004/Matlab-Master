%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Automatic Conversion of Hamamatsu TIFF images into a tiff stack%%%%%%%%%%
%%AJN 5-30-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all

%Select the first image in your future stack. The Hamamatsu outputs files
%in a sequential pattern that this program takes advantage of.
[sname spath] = uigetfile('.Tif' , ' Choose the first image in your series.');
stack_name = sname(1:numel(sname)-6);
fpath = strcat(spath,sname,'.tif');
destpath = 'J:\Data\5-31-13 CMOS Dendra2\';  %destination path

cd(spath);

l=1;


k=1;

%Finishes the image array
while l == 1
    sname1 = strcat(stack_name,'_', num2str(k));             %Picks next file in series
    fpath1 = strcat(spath , sname1, '.tif');   %Builds full path + file
    if exist(fpath1 , 'file') == 2                  %Checks for existance of file
        im1(:,:,k) = imread(fpath1, 'Tiff');        %Writes file to array if exists
        k=k+1;                  
    else                                            %Assumes we have reached the end of the file list
        l=0;
    end
    
end


%Writes the array as a single tiff
outfname = strcat(destpath,stack_name);

save(outfname, 'im1');