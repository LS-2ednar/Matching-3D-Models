%% Setup reading files and creating point clouds
clear;clc;close all;

% reading stl files
stlData = stlread('Mand-left-cut.stl');
mand = stlData.Points;
stlData1 = stlread('Pelvis-left-cut.stl');
pelvis = stlData1.Points;

% Initial Overview 3D plots of both stl objects sepeartle
figure
mand_fig = plot3(mand(:,1),mand(:,2),mand(:,3),'.');
title('mand')

figure
pelvis_fig = plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
title('pelvis')

%plot of both in one figure
figure
title('mixed')
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
hold on
pelvis_fig = plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');