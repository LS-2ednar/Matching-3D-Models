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
mand = move(mand,pelvis);
% 3d plot of both parts intially
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
hold on
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
title('mixed')

%% Translating and Rotating the mand

index = model_iteration(mand,pelvis,3,1); %mand, pelvis, number of subsets, size scaling of second point cloud 


for i = 1:(size(index,2)-1)
    fprintf('Subset %d\n' ,i);
    %select part from pelvis which shall be matched by mand
    figure
    if i == 1
        plot3(pelvis(index(i):index(i+1),1),pelvis(index(i):index(i+1),2),pelvis(index(i):index(i+1),3),'.b'); %mainpart
        hold on
        plot3(pelvis(index(i+1):end,1),pelvis(index(i+1):end,2),pelvis(index(i+1):end,3),'.','Color','#AAAAAA'); %mainpart
        new_pelvis = pelvis(index(i):index(i+1),:);
        title(['Subset: ',num2str(i)])
%     elseif i == size(index,2)
%         plot3(pelvis(index(i-1):index(i),1),pelvis(index(i-1):index(i),2),pelvis(index(i-1):index(i),3),'.b'); %mainpart
%         hold on
%         plot3(pelvis(1:index(i),1),pelvis(1:index(i),2),pelvis(1:index(i),3),'.','Color','#AAAAAA')%bevor
%         new_pelvis = pelvis(index(i-1):index(i),:)
%         title('case 2')
    else
        plot3(pelvis(index(i):index(i+1),1),pelvis(index(i):index(i+1),2),pelvis(index(i):index(i+1),3),'.b'); %main part
        hold on
        plot3(pelvis(1:index(i),1),pelvis(1:index(i),2),pelvis(1:index(i),3),'.','Color','#AAAAAA')%bevor
        plot3(pelvis(index(i+1):end,1),pelvis(index(i+1):end,2),pelvis(index(i+1):end,3),'.','Color','#AAAAAA')%after
        new_pelvis = pelvis(index(i):index(i+1),:);
        title(['Subset: ',num2str(i)])
        
    end
    %determine new starting point for mand
    new_mand = move(mand, new_pelvis);
    plot3(new_mand(:,1),new_mand(:,2),new_mand(:,3),'.m');
    drawnow
    
    %Optimization procedure
    
    %implementation of Simulated Annealing algorithm
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
    rejected = [0, 0, 0, 0, 0, 0];
    
    % Use hausdorff distance of the inital positions as initial best value
	distance_best = directed_hausdorff_distance(new_mand, new_pelvis);
    
    % Calculate boundaries for the solution space
    x_max = max(new_pelvis(:,1));
    x_min = min(new_pelvis(:,1));
    y_max = max(new_pelvis(:,2));
    y_min = min(new_pelvis(:,2));
    z_max = max(new_pelvis(:,3));
    z_min = min(new_pelvis(:,3));
    
    % Set starting temperature for the outer loop, the max stepsize and the max
    % rotation
    startT = 10;
    maxStep = 10;
    maxRotation = 1;
    
    fprintf('Start Fitting of Subset %d of %d\n', i, (size(index,2)-1)); 
   
    while distance_best > 10^(-2)
        for T=startT:-1:1
            fprintf([num2str(T),' more Subruns'])
            for v=1:10
                fprintf(['Subrun: ',num2str(v),' of 10\n'])
                % randomly update parameters for rotation
                parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation*T/startT;
                parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation*T/startT;
                parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation*T/startT;

                % randomly update parameters for translation
                parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
                parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
                parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;

                % transform the mand matrix
                mand_current = transformation(parameters_current, new_mand);

                % check if parameters were already rejected 
                tf = ismember(parameters_current, rejected, 'rows');

                % update the parameters as long as we are not in the solution space
                % or are already in the rejected parameters
                while (max(mand_current(:,1)) >= x_max || min(mand_current(:,1)) <= x_min || ...
                       max(mand_current(:,2)) >= y_max || min(mand_current(:,2)) <= y_min || ...
                       max(mand_current(:,3)) >= z_max || min(mand_current(:,3)) <= z_min || ...
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
                distance_current = hausdorff_distance(mand_current, new_pelvis);
                difference = distance_current - distance_best;

                % if the new distance is smaller than the last distance accept the
                % solution
                if difference < 0
                    parameters_best = parameters_current;
                    distance_best = distance_current;

                % else if the new distance is not smaller than the last distance
                % accept the solution with a random probability
                elseif (exp((-difference*30)/T) > rand)
                    parameters_best = parameters_current;
                    distance_best = distance_current;
                else
                    rejected = [rejected; parameters_current];
                end

            end
            plot3(mand_current(:,1),mand_current(:,2),mand_current(:,3),'.')
            distance_best
            fprintf(['best distance: 'num2str(distance_best)])
            drawnow
        end
    end

    distance_best
end