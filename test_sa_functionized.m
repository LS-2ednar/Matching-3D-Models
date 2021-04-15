%---TEST SECTION-----------------------------------------------------------
clear; clc; close all;
fprintf('Reading files\n')
stlData = stlread('Mand-left-cut.stl');
mand = stlData.Points;
stlData1 = stlread('Pelvis-left-cut.stl');
pelvis = stlData1.Points;

fprintf('Start Simulated Anneling\n')
[new_fig, bd] = simulatedAnneling(mand,pelvis,[5,5],[50,10],[1,0.1],[1,0.1])

fprintf('Ploting')
figure
pelvis_fig = plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
hold on
plot3(new_fig(:,1),new_fig(:,2),new_fig(:,3),'.g')
title('simulated anneling')
%---FUNCTION---------------------------------------------------------------
function [end_mand, end_distances] = simulatedAnneling(mand,pelvis,number_of_restarts,startT,maxStep,maxRotation)
%Two Step Simulated Anneling with an inital rough approximation flowed by a
%fine tuning of sveral good inital results
%Parameters are: mand = smaller point colud, pelvis = biger point cloud,
%number_of_restarts = number of generated results array[1x2], startTemp = Inital
%Temperature array[1x2], maxStep = maximal Step size array [1x2], maxRotation = maximal number of
%roations array[1x2]

%% Simulated Annealing algorithm
% Initialize paramters Rotation matrix to unit matrix and translation vector 
% to zero vector

alpha = 0;
beta = 0;
gamma = 0;
xt = 0;
yt = 0;
zt = 0;

parameters_best = [alpha, beta, gamma, xt, yt, zt];
parameters_current = parameters_best;

% Create matrix to remember rejected solutions to get a shorter running
% time, since the calculation of the (modified) hausdorff distance is 
% relatively time consuming and add the inital parameters 
% rejected = parameters_current;

%%
% Use hausdorff distance of the inital positions as initial best value
fprintf('Initial distance calculation: ')
tic
distance_best = directed_averaged_hausdorff_distance(mand, pelvis, 10);
toc
%%
% Calculate boundaries for the solution space
x_max = max(pelvis(:,1));
x_min = min(pelvis(:,1));
y_max = max(pelvis(:,2));
y_min = min(pelvis(:,2));
z_max = max(pelvis(:,3));
z_min = min(pelvis(:,3));

% Set starting temperature for the outer loop, the max stepsize and the max
% rotation
% number_of_restarts = 5;
% startT = 50;
% maxStep = 1;
% maxRotation = 1;
best_positions = [];
best_distances = [];
fprintf('Find inital Solutions\n')
%initial solution
for restarts=1:number_of_restarts(1)
    fprintf(['Working on ',num2str(restarts),' restart from ',num2str(number_of_restarts(1)),'.\n'])
    for T=startT(1):-1:1

        for v=1:5
            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation(1)*T/startT(1);
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation(1)*T/startT(1);
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation(1)*T/startT(1);

            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep(1)*T/startT(1);
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep(1)*T/startT(1);
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep(1)*T/startT(1);

            % transform the mand matrix
            mand_current = transformation(parameters_current, mand);

            % update the parameters as long as we are not in the solution space
            % or are already in the rejected parameters
            while (max(mand_current(:,1)) > x_max+10 || min(mand_current(:,1)) < x_min-10 || ...
                   max(mand_current(:,2)) > y_max+10 || min(mand_current(:,2)) < y_min-10 || ...
                   max(mand_current(:,3)) > z_max+10 || min(mand_current(:,3)) < z_min-10)
               
            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation(1)*T/startT(1);
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation(1)*T/startT(1);
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation(1)*T/startT(1);

            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep(1)*T/startT(1);
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep(1)*T/startT(1);
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep(1)*T/startT(1);

            % transform the mand matrix
            mand_current = transformation(parameters_current, mand);

            % check if parameters were already rejected 
            %tf = ismember(parameters_current, rejected, 'rows');
            end

            % calculated the (modified) hausdorff distance for the transformed
            % mand matrix 
            distance_current = directed_averaged_hausdorff_distance(mand_current, pelvis,10);
            difference = distance_current - distance_best;

            % if the new distance is smaller than the last distance accept the
            % solution
            if difference < 0
                parameters_best = parameters_current;
                distance_best = distance_current;

            % else if the new distance is not smaller than the last distance
            % accept the solution with a random probability
            elseif (exp((-difference*300)/T) > rand)
                parameters_best = parameters_current;
                distance_best = distance_current;
            fprintf('Works')
            end
        fprintf(['Temperature: ',num2str(T)])
        end
        
    end
    best_positions = [best_positions; parameters_best];
    best_distances = [best_distances; distance_best];
    
end

%angeben der besten distanzen ()

%%
%fine tuning
fprintf('Fine tuning for final Solutions\n')
end_mand = {};
end_distances = [];

for i=1:number_of_restarts(2)
    mand_new = transformation(best_positions(i,:), mand);

%     figure,
%     plot3(mand_new(:,1),mand_new(:,2),mand_new(:,3),'.')
%     hold on
%     plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
%     title('Fine Tuning')

    % Set starting temperature for the outer loop, the max stepsize and the max
    % rotation
    alpha = 0;
    beta = 0;
    gamma = 0;
    xt = 0;
    yt = 0;
    zt = 0;

    parameters_best = [alpha, beta, gamma, xt, yt, zt];
    parameters_current = parameters_best;

%     startT = 10;
%     maxStep = 0.1;
%     maxRotation = 0.1;
    fprintf(['Working on ',num2str(restarts),' restart from ',num2str(number_of_restarts(2)),'.\n'])
    for T=startT(2):-1:1

        for v=1:5
            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation(2)*T/startT(2);
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation(2)*T/startT(2);
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation(2)*T/startT(2);
            
            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep(2)*T/startT(2);
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep(2)*T/startT(2);
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep(2)*T/startT(2);
            
            % transform the mand matrix
            mand_current = transformation(parameters_current, mand_new);

            % update the parameters as long as we are not in the solution space
            % or are already in the rejected parameters
            while (max(mand_current(:,1)) > x_max || min(mand_current(:,1)) < x_min || ...
                   max(mand_current(:,2)) > y_max || min(mand_current(:,2)) < y_min || ...
                   max(mand_current(:,3)) > z_max || min(mand_current(:,3)) < z_min)

            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation(2)*T/startT(2);
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation(2)*T/startT(2);
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation(2)*T/startT(2);
            
            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep(2)*T/startT(2);
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep(2)*T/startT(2);
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep(2)*T/startT(2);

            % transform the mand matrix
            mand_current = transformation(parameters_current, mand_new);

            end

            % calculated the (modified) hausdorff distance for the transformed
            % mand matrix 
            distance_current = directed_averaged_hausdorff_distance(mand_current, pelvis,5);
            difference = distance_current - distance_best;

            % if the new distance is smaller than the last distance accept the
            % solution
            if difference < 0
                parameters_best = parameters_current;
                distance_best = distance_current;

            % else if the new distance is not smaller than the last distance
            % accept the solution with a random probability
            elseif (exp((-difference*300)/T) > rand)
                parameters_best = parameters_current;
                distance_best = distance_current;

            end
        end
%         plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
%         drawnow
        fprintf(['Temperature: ',num2str(T)])
    end
    end_mand{i} = transformation(parameters_best, mand_new);
    end_distances = [end_distances; distance_best];
    
end

end






% 
% %--------------------NICHT HINZUFüGEN SOLLTE SCHON VORHANDEN SEIN----------
% function [X_new] = transformation(parameters, X)
% %Transform a point cloud in 3D
% alpha = parameters(1);
% beta = parameters(2);
% gamma = parameters(3);
% xt = parameters(4);
% yt = parameters(5);
% zt = parameters(6);
% r11 = cos(alpha)*cos(beta);
% r12 = cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma);
% r13 = cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
% r14 = xt;
% 
% r21 = sin(alpha)*cos(beta);
% r22 = sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma);
% r23 = sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma);
% r24 = yt;
% 
% r31 = -sin(beta);
% r32 = cos(beta)*sin(gamma);
% r33 = cos(beta)*cos(gamma);
% r34 = zt;
% 
% r41 = 0;
% r42 = 0;
% r43 = 0;
% r44 = 1;
% 
% T = [r11, r12, r13, r14;...
%     r21, r22, r23, r24;...
%     r31, r32, r33, r34;...
%     r41, r42, r43, r44];
% 
% extension = ones(length(X),1);
% X_ext = horzcat(X, extension);
% X_new_ext = T*X_ext';
% X_new = X_new_ext(1:3,:)';
% end
