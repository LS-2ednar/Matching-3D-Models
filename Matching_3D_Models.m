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

%% Translating and Rotating the mand

alpha = 0;
beta = 0;
gamma = 0;
xt = 0;
yt = 0;
zt = 0;

mand_new = mand;

function [mand_new, distance] = ObjectiveFunction(alpha, beta, gamma, xt, yt, zt, mand_new, pelvis)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

r11 = cos(alpha)*cos(beta);
r12 = cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma);
r13 = cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
r14 = xt;

r21 = sin(alpha)*cos(beta);
r22 = sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma);
r23 = sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma);
r24 = yt;

r31 = -sin(beta);
r32 = cos(beta)*sin(gamma);
r33 = cos(beta)*cos(gamma);
r34 = zt;

r41 = 0;
r42 = 0;
r43 = 0;
r44 = 1;

T = [r11, r12, r13, r14;...
    r21, r22, r23, r24;...
    r31, r32, r33, r34;...
    r41, r42, r43, r44];

extension = ones(length(mand),1);
mand_ext = horzcat(mand, extension);
mand_new_ext = T*mand_ext';
mand_new = mand_new_ext(1:3,:)';

distance = hausdorff(mand_hausdorff', pelvis);
end
