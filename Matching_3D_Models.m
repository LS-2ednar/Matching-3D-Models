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
%updating mand position to be closer to the pelvis
mand = move(mand,300,80,300);
% 3d plot of both parts intially
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
hold on
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
% xlabel('x')
% ylabel('y')
% zlabel('z')
title('mixed')

%% Translating and Rotating the mand

alpha = 0;
beta = 0;
gamma = 0;
xt = 0;
yt = 0;
zt = 0;

mand_new = mand;

% Lukas Gebastell
[mand_new, distance] = ObjectiveFunction(alpha,beta,gamma,xt,yt,zt,mand_new,pelvis)
figure
plot3(mand_new(:,1),mand_new(:,2),mand_new(:,3),'.')
hold on
pelvis_fig = plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
title('mixed with new mand')

for i = 1:100
    alpha = alpha+1;
    beta = beta+1;
    gamma = gamma+1;
    xt = xt+1;
    yt = yt+1;
    zt = zt+1;
    [mand_new, distance] = ObjectiveFunction(alpha,beta,gamma,xt,yt,zt,mand_new,pelvis);
    plot3(mand_new(:,1),mand_new(:,2),mand_new(:,3),'.')
    
end

mand = mand_new; % remove later since this is just here for fun :-D
% Optimization


% Functions
% ObjectiveFunction
function [mand_new, distance] = ObjectiveFunction(alpha, beta, gamma, xt, yt, zt, mand, pelvis)
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

distance = hausdorff(mand_new', pelvis);
end

%Utility Function to move initial plot closer to the actual problem
function stl = move(stl,x,y,z)
stl(:,1) = stl(:,1) + x;
stl(:,2) = stl(:,2) + y;
stl(:,3) = stl(:,3) + z;
new_stl = stl;
end