%% Setup reading files and creating point clouds
clear;clc;close all;

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
plot3(mand(:,1),mand(:,2),mand(:,3),'m.')
hold on
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
% xlabel('x')
% ylabel('y')
% zlabel('z')
title('Trial')
drawnow
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

%% Use hausdorff distance of the inital positions as initial best value
fprintf('Initial Alignment: ')
tic
distance_best = directed_averaged_hausdorff_distance(mand, pelvis,10);
toc
%% Calculate boundaries for the solution space
x_max = max(pelvis(:,1));
x_min = min(pelvis(:,1));
y_max = max(pelvis(:,2));
y_min = min(pelvis(:,2));
z_max = max(pelvis(:,3));
z_min = min(pelvis(:,3));

% Set starting temperature for the outer loop, the max stepsize and the max
% rotation
startT = 100;
maxStep = 5;
maxRotation = 1;
store_attempts = 5;

%initialization of later potential attempts
potential_attempts = rand(store_attempts,7)*1000; % 1 = distance, 2-7 = parameters of transitionmatrix

while distance_best > 5*10^(-1)
    tic
    fprintf('Start Simulated Anneling: Quick Part\n')
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
            
            % check if parameters were already rejected 
            tf = ismember(parameters_current, rejected, 'rows');
            
            % update the parameters as long as we are not in the solution space
            % or are already in the rejected parameters
            while (max(mand_current(:,1)) > x_max+5 || min(mand_current(:,1)) < x_min-5 || ... 
                   max(mand_current(:,2)) > y_max+5 || min(mand_current(:,2)) < y_min-5 || ...
                   max(mand_current(:,3)) > z_max+5 || min(mand_current(:,3)) < z_min-5 || ...
                   tf)                                                                                  

            % record rejected parameters
            rejected = [rejected; parameters_current];

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
            tf = ismember(parameters_current, rejected, 'rows');
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
            elseif (exp((-difference*50)/T) > rand)                                                         
                p = exp((-difference*50)/T)                                                                 
                parameters_best = parameters_current;
                distance_best = distance_current;
            end
            rejected = [rejected; parameters_current];
            
        end

        plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
        T
        distance_best
        if distance_current < min(potential_attempts(:,1))
            %maximum anpasse in matrix
            idx = potential_attempts == max(potential_attempts(:,1)); % get index of max distance and replace row afterwards
            index = find(idx);
            potential_attempts(index,1) = distance_current;
            potential_attempts(index,2) = parameters_current(1);
            potential_attempts(index,3) = parameters_current(2);
            potential_attempts(index,4) = parameters_current(3);
            potential_attempts(index,5) = parameters_current(4);
            potential_attempts(index,6) = parameters_current(5);
            potential_attempts(index,7) = parameters_current(6);
        end
%         parameters_best
%         parameters_current
        drawnow
        
    end
    toc
%     fprintf(num2str(store_attempts)+' best where estimated\n\n')
    
end
%% Ploting n best Results after initial plot

fprintf('Ploting Inital Potential Results\n\n')
figure
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
hold on
title('Fit after initial rough optimization')
for i = 1:size(potential_attempts,1)
    if potential_attempts(i,1) < 1
        mand_current = transformation(potential_attempts(i,2:7), mand);
        plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
        drawnow
    end
end

%% refinement of the preselected attempts

figure
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
hold on
title('Refinement')
for i = 1:size(potential_attempts,1)
    parameters_current = potential_attempts(i,2:7)
    if potential_attempts(i,1) < 1
        while distance_best > 5*10^(-3)
            tic
            fprintf('Start Simulated Anneling: Quick Part\n')
            for T=startT:-1:1

                for v=1:5
                    % randomly update parameters for rotation
                    parameters_current(1) = potential_attempts(i,2) + (rand-0.5)*2*maxRotation*T/startT;
                    parameters_current(2) = potential_attempts(i,3) + (rand-0.5)*2*maxRotation*T/startT;
                    parameters_current(3) = potential_attempts(i,4) + (rand-0.5)*2*maxRotation*T/startT;

        %             % randomly update parameters for translation
        %             parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
        %             parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
        %             parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;

                    % transform the mand matrix
                    mand_current = transformation(parameters_current, mand);

                    % check if parameters were already rejected 
                    tf = ismember(parameters_current, rejected, 'rows');

                    % update the parameters as long as we are not in the solution space
                    % or are already in the rejected parameters
                    while (max(mand_current(:,1)) > x_max+5 || min(mand_current(:,1)) < x_min-5 || ... 
                           max(mand_current(:,2)) > y_max+5 || min(mand_current(:,2)) < y_min-5 || ...
                           max(mand_current(:,3)) > z_max+5 || min(mand_current(:,3)) < z_min-5 || ...
                           tf)                                                                                  %maxStep einsetzen?

                    % record rejected parameters
                    rejected = [rejected; parameters_current];

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
                    tf = ismember(parameters_current, rejected, 'rows');
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
                    elseif (exp((-difference*50)/T) > rand)                                                         %tmax anstelle von 50?
                        p = exp((-difference*50)/T)                                                                 %tmax anstelle von 50?
                        parameters_best = parameters_current;
                        distance_best = distance_current;
                    end
                    rejected = [rejected; parameters_current];

                end

                plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
                T
                distance_best
        %         parameters_best
        %         parameters_current
                drawnow

            end
            toc
        %     fprintf(num2str(store_attempts)+' best where estimated\n\n')

        end
    end
end

%plot best solution
figure
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
hold on
title('Best Solution')
plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
