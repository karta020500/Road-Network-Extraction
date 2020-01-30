close all;
clear all;

[bundleImage,R] = geotiffread('Road_calss.tif');

[H,theta,rho] = hough(bundleImage);
P = houghpeaks(H,20,'threshold',0);
lines = houghlines(bundleImage,theta,rho,P,'FillGap',100,'MinLength',50);
H = lines;


