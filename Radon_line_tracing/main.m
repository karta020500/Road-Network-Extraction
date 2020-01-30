close all;
clear all;
nFilesToLookAt=0; % 0 for all
folder_name = uigetdir;
% Get list of all *.png files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
imagefiles = dir(fullfile(folder_name,'/*.tif'));
nfiles = length(imagefiles);
if nFilesToLookAt==0
    nFilesToLookAt=nfiles;
end

AllWidths=[];
for i=1:round(nfiles/nFilesToLookAt):nfiles
    currentfilename = fullfile(folder_name,imagefiles(i).name);
    [w, long_seg] = linewidthofimages(currentfilename);
    AllWidths=[AllWidths w];
end 

figure;
plot(AllWidths);

disp(sprintf('mean width:%0.2f',mean(AllWidths)));
disp(sprintf('median width:%0.2f',median(AllWidths)));