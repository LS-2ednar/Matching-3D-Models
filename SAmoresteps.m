%% Setup reading files and creating point clouds
clear;clc;close all;
tic
% reading stl files
stlData = stlread('Mand-left-cut.stl');
mand = stlData.Points;
stlData1 = stlread('Pelvis-left-cut.stl');
pelvis = stlData1.Points;

% Initial overview 3D plots of both stl objects 
figure
mand_fig = plot3(mand(:,1),mand(:,2),mand(:,3),'.');
title('mand')
figure
pelvis_fig = plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
title('pelvis')

% updating mand position, the mand point cloud is moved to the center of
% gravity of the pelvis point cloud
mand = move(mand,pelvis);

%plot of both in one figure
figure
% 3d plot of both parts intially
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
hold on
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
% xlabel('x')
% ylabel('y')
% zlabel('z')
title('mixed')

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
rejected = parameters_current;

%%
% Use hausdorff distance of the inital positions as initial best value
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
number_of_restarts = 5;
startT = 50;
maxStep = 1;
maxRotation = 1;
best_positions = [];
best_distances = [];

for restarts=1:number_of_restarts
    for T=startT:-1:1

        for v=1:5
            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation*T/startT;

            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;

            % transform the mand matrix
            mand_current = transformation(parameters_current, mand);

            % update the parameters as long as we are not in the solution space
            % or are already in the rejected parameters
            while (max(mand_current(:,1)) > x_max+10 || min(mand_current(:,1)) < x_min-10 || ...
                   max(mand_current(:,2)) > y_max+10 || min(mand_current(:,2)) < y_min-10 || ...
                   max(mand_current(:,3)) > z_max+10 || min(mand_current(:,3)) < z_min-10)
               
            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation*T/startT;

            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;

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

            end

        end
    end
    best_positions = [best_positions; parameters_best];
    best_distances = [best_distances; distance_best];
end


distance_best



%%
%plot of both in one figure
figure,
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
hold on
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
title('Results from multiple restarts')

for i=1:number_of_restarts

    mand_current = transformation(best_positions(i,:), mand);
    plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')

end
hold off
toc

%%
%fine tuning
end_mand = {};
end_distances = [];

for i=1:number_of_restarts
    mand_new = transformation(best_positions(i,:), mand);

    figure,
    plot3(mand_new(:,1),mand_new(:,2),mand_new(:,3),'.')
    hold on
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
    title('Fine Tuning')

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

    startT = 10;
    maxStep = 0.1;
    maxRotation = 0.1;

    for T=startT:-1:1

        for v=1:5
            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation*T/startT;
            
            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;
            
            % transform the mand matrix
            mand_current = transformation(parameters_current, mand_new);

            % update the parameters as long as we are not in the solution space
            % or are already in the rejected parameters
            while (max(mand_current(:,1)) > x_max || min(mand_current(:,1)) < x_min || ...
                   max(mand_current(:,2)) > y_max || min(mand_current(:,2)) < y_min || ...
                   max(mand_current(:,3)) > z_max || min(mand_current(:,3)) < z_min)

            % randomly update parameters for rotation
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation*T/startT;
            
            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;

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
        plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
        drawnow
    end
    end_mand{i} = transformation(parameters_best, mand_new);
    end_distances = [end_distances; distance_best];
end


%%
%%
%plot of both in one figure
figure,
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
hold on
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
title('Results from multiple restarts')

for i=1:number_of_restarts

    mand_current = end_mand{1,i};
    plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')

end
hold off
toc
