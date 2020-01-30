close all
clear
clc

%% get demo image
% bw = imread('demo.png');
% bw = im2bw(bw);
% bw = imdilate(bw, strel('disk', 1));

[bundleImage,R] = geotiffread('Road_calss.tif');
R2 = worldfileread('Road_calss.tfw', 'planar', size(bundleImage));

bundleImage(bundleImage > 3) = 0;
bundleImage(bundleImage == 3) = 1;
bundleImage=double(bundleImage);

%pre-processing------------------------------------------------------------

%opening
se = strel('square',2);
opened = imopen(bundleImage,se);

%closeing
se2 = strel('square',5);
closed = imclose(opened,se2);
%closed = imclose(closed,se2);
%connectivity_area
Pixel_connectivity = bwareaopen(closed, 10);


%hole_fill
filled = imfill(Pixel_connectivity,'holes');
holes = filled & ~Pixel_connectivity;
bigholes = bwareaopen(holes, 200);
smallholes = holes & ~bigholes;
filled = Pixel_connectivity | smallholes;

skel = bwmorph(filled,'skel',Inf);

skel = bwmorph(skel, 'spur', 9);

skel = imdilate(skel,ones(5,5));

block_image = filled;
%% set up parameters for line detection

p = struct();
p.theta_step = 0.5; % resolution for theta values (selectivity of line angles)
p.peak_threshold = 0.3; % threshold for peak detection (0.3 - 30% of max intensity)
p.line_fillgap = 60; % if two adjucent lines have distance smaller - merge them
p.line_minlength = 10; % if line is smaller than this, do not count it
p.split_bw = 0; % split image into segments or not
p.show_bw_plots = 1; % show plot with found lines for first 20 segments (for debugging)
p.angletol = 0; % if two lines cross each other with angle smaller than this, merge them

%% detect lines
bwlines = getlinesforbw(block_image, p);
disp(bwlines);

%% show the detected lines
figure
showfoundlines(block_image, bwlines)