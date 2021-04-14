%% -------- SETUP: Load Data & Create and Visualize Point Clouds ---------
clear;clc;close all;
tic
% reading stl files
stlData = stlread('Mand-left-cut.stl');
mand = stlData.Points;
stlData1 = stlread('Pelvis-left-cut.stl');
pelvis = stlData1.Points;

% Initial overview 3D plots of both stl objects 
figure(1),
plot3(mand(:,1),mand(:,2),mand(:,3),'.');
title('Initial Point Cloud: Mandible')

figure(2),
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
title('Inital Point Cloud: Pelvis')

% updating mand position, the mand point cloud is moved to the center of
% gravity of the pelvis point cloud
mand = move(mand,pelvis);

%plot of both in one figure
figure(3),
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
hold on
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
title('The inital position of the mandible and pelvis')
hold off

%% ------------------ OPTIMIZATION: SIMULATED ANNEALING ------------------

%% -------------------- OPTIMIZATION: PARTICLE SWARM ---------------------
% Quick Alignment Parameters
MinRot_Q = 0;            % minimum value for rotation parameters 
MaxRot_Q = 2*pi;         % maximum value for rotation parameters
MinTrans_Q = -100;       % minimum value for translation parameters
MaxTrans_Q = 100;        % maximum value for translation parameters
iter_Q = 50;             % number of iterations
nPop_Q = 20;             % Population Size (Swarm Size)
wdamp_Q = 0.98;          % damping coefficient
c1_Q = 1.5;              % personal acceleration coefficent
c2_Q = 1.5;              % social acceleration coefficient
stepsize_Q = 10;         % stepsize for the hausdorff distance

% Fine Alignment Parameters
MinRot_F = 0;            % minimum value for rotation parameters 
MaxRot_F = 2*pi;         % maximum value for rotation parameters
MinTrans_F = -10;        % minimum value for translation parameters
MaxTrans_F = 10;         % maximum value for translation parameters
iter_F = 50;             % number of iterations
nPop_F = 5;              % Population Size (Swarm Size)
wdamp_F = 0.98;          % damping coefficient
c1_F = 1.1;              % personal acceleration coefficent
c2_F = 1.1;              % social acceleration coefficient
stepsize_F = 5;          % stepsize for the hausdorff distance

restarts = 5;            % number of restarts

mandible.Points = mand;
mandible.Distance = inf;

% save the best mand and its distance of the quick alignment for each restart 
MandQuick = repmat(mandible, restarts,1); 
% save the best mand and its distance of the fine alignment for each restart 
MandFine = repmat(mandible, restarts,1);

for i=1:5
    
    % quick alignment
    [mand_quick, distances_quick] = ParticleSwarmOpti(MinRot_Q, MaxRot_Q, MinTrans_Q, ...
        MaxTrans_Q, iter_Q, nPop_Q, wdamp_Q, c1_Q, c2_Q, mand, pelvis, stepsize_Q);
    
    MandQuick(i).Points = mand_quick;
    MandQuick(i).Distance = min(distances_quick);
    
    figure(4),
    plot(distances_quick, 'Linewidth', 1);
    xlabel('Iteration')
    ylabel('Best Distance')
    title('Best Distance per Iteration for Quick Alignment')
    hold on
    drawnow
    
    figure(5),
    hold on
    plot3(mand_quick(:,1),mand_quick(:,2),mand_quick(:,3),'.')
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
    title('Best Transformations for Quick Alignment')
    drawnow
    
    % fine alignment
    [mand_fine, distances_fine] = ParticleSwarmOpti(MinRot_F, MaxRot_F, MinTrans_F, ...
        MaxTrans_F, iter_F, nPop_F, wdamp_F, c1_F, c2_F, mand_quick, pelvis, stepsize_F);
    
    MandFine(i).Points = mand_fine;
    MandFine(i).Distance = min(distances_fine);
    
    figure(6),
    plot(distances_fine, 'Linewidth', 1);
    xlabel('Iteration')
    ylabel('Best Distance')
    title('Best Distance per Iteration for Fine Alignment')
    hold on
    drawnow
    
    figure(7),
    hold on
    plot3(mand_fine(:,1),mand_fine(:,2),mand_fine(:,3),'.')
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
    title('Best Transformations for Fine Alignment')
    drawnow

end
hold off
toc
%% ------------------------------ FUNCTIONS ------------------------------

function [X] = move(X,Y)
% moves the first point set X to the center of gravity of the second point set
% INPUT
% X: point set
% Y: point set
% OUTPUT:
% X: moved point set X

dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    fprintf('All points must have the same dimension')
end

% calculate the center of gravity for both point sets
center_of_gravity_X = zeros(1,dimX(2));
center_of_gravity_Y = zeros(1,dimY(2));

for i=1:dimX(2)
    x = sum(X(:,i))/dimX(1);
    y = sum(Y(:,i))/dimY(1);
    center_of_gravity_X(i) = x;
    center_of_gravity_Y(i) = y;
end

% calculate the distance between the two centers of gravity
dist_vector = abs(center_of_gravity_X - center_of_gravity_Y);

% move the first point set
for i=1:dimX(2)
    X(:,i) = X(:,i) + dist_vector(i);
end

end

function [X_new] = transformation(parameters, X)
% Transforms a 3D point set
% INPUT:
% parameters: a 1x6 vector with entries corresponding to the rotation and
%             translation values
% X:          a 3D point set
% OUTPUT:
% X_new:      the transformed 3D point set

dimX = size(X);
dim_p = size(parameters);

if dimX(2) ~= 3
    fprintf('Input has to be a 3D point set')
end

if dim_p(2) ~= 6
    fprint('Function needs one parameter vector with 6 entries')
end

alpha = parameters(1);  % rotation around the x-axis
beta = parameters(2);   % rotation around the y-axis
gamma = parameters(3);  % rotation around the z-axis
xt = parameters(4);     % translation along the x-axis 
yt = parameters(5);     % translation along the y-axis
zt = parameters(6);     % translation along the z-axis

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

% rigid transformation matrix
T = [r11, r12, r13, r14;...
    r21, r22, r23, r24;...
    r31, r32, r33, r34;...
    r41, r42, r43, r44];

extension = ones(length(X),1);
X_ext = horzcat(X, extension);
X_new_ext = T*X_ext';
X_new = X_new_ext(1:3,:)';
end

function [dahd] = directed_averaged_hausdorff_distance(X,Y, step)
% Calculates the directed hausdorff distance from X to Y 
% If the directed averaged Hausdorff distance is zero all points of X lie on a
% point in Y, this can only happen if X <= Y.

% INPUT:
% X: a 3D point set
% Y: a 3D point set 
% OUTPUT:
% dahd: the directed averaged hausdorff distance from X to Y

dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    fprintf('All points must have the same dimension')
end

if dimX(1) > dimY(1)
    fprint('The first input point set is larger than the second one - hausdorff distance cannot reach zero')
end

%calculate the directed distance from X to Y
dXY_all = zeros(1, dimX(1));
for i=1:step:dimX(1)
    shortestdist = norm(X(i,:) - Y(1,:));
    for j=2:step:dimY(1)
        dist = norm(X(i,:) - Y(j,:));
        if dist < shortestdist
            shortestdist = dist;
        end
    end
    dXY_all(i) = shortestdist;
end
% average the distance over all points
dahd = sum(dXY_all)/(dimX(1)/step);
end


function [ObjectMoveNew, BestDistances] = ParticleSwarmOpti(MinRotation, MaxRotation,...
    MinTranslation, MaxTranslation, iter, nPopulation, wdamp, personalacceleration, socialacceleration, ObjectMove, ObjectFixed, stepsize)
% Function for a particle swarm optimization for registration of two 3D 
% point clouds - find the best transformation of ObjectMove which minimized
% the directed average hausdorff distance to ObjectFixed

% INPUT:
% MinRotation:          minimum value for rotation parameters
% MaxRotation:          maximum value for rotation parameters
% MinTranslation:       minimum value for translation parameters
% MaxTranslation:       maximum value for translation parameters
% iter:                 number of iterations
% nPopulation:          number of particles in the population
% wdamp:                damping parameter for the velocity
% personalacceleration: personal acceleration coefficient
% socicalacceleration:  social acceleration coeffficient
% ObjectMove:           object which is moved though the soultion space
% ObjectFixed:          object which we want ObjectMove to align to
% stepsize:             gives the length of the interval between points to 
%                       evaluate in the directed averaged hausdorff distance
% OUTPUT:
% ObjectMoveNew:        the moved 3D point cloud
% BestDistance:         the distance of the new point cloud to the fixed 
%                       point cloud

% The Objective Function we want to minimize
% x: vector of transformation parameters
ObjFunc = @(x) directed_averaged_hausdorff_distance(transformation(x, ObjectMove), ObjectFixed, stepsize);

nParameters = 6;    % number of transformation parameters
w = 1;              % Inertia coefficient

% initialize the transformation parameters to zero
% and Velocity and Distance as empty arrays

inital_particle.Transformation = zeros(1,nParameters);
inital_particle.Velocity = [];
inital_particle.Distance = [];
inital_particle.Best.Transformation = [];
inital_particle.Best.Distance = [];

particle = repmat(inital_particle, nPopulation,1); % initialize the whole population

GlobalBest.Distance = inf;                       % initialize the worst case


for i=1:nPopulation
    
    % set random values for the rotation parameters between MinRotation and
    % MaxRotation
    particle(i).Transformation(:,1:3) = unifrnd(MinRotation, MaxRotation, [1 3]);
    % set random values for the translation parameters between
    % MinTranslation and MaxTranslation
    particle(i).Transformation(:,4:6) = unifrnd(MinTranslation, MaxTranslation, [1 3]);
    
    % initialize the velocity to zero
    particle(i).Velocity = zeros([1 nParameters]);
    
    particle(i).Distance = ObjFunc(particle(i).Transformation);
    
    particle(i).Best.Transformation = particle(i).Transformation;
    particle(i).Best.Distance = particle(i).Distance;
    
    if particle(i).Best.Distance < GlobalBest.Distance
        GlobalBest = particle(i).Best;
    end

end

% hold best global distance value for each iteration
BestDistances = zeros(iter, 1);

for it=1:iter
    
    for i=1:nPopulation
        % update velocity for each particle 
        particle(i).Velocity = w*particle(i).Velocity + ...
            personalacceleration*rand([1 nParameters]).*(particle(i).Best.Transformation - particle(i).Transformation) + ...
            socialacceleration*rand([1 nParameters]).*(GlobalBest.Transformation - particle(i).Transformation);
        
        % update transformation parameters with the new velocity 
        particle(i).Transformation = particle(i).Transformation + particle(i).Velocity;
        
        % evaluate the ObjFunc
        particle(i).Distance = ObjFunc(particle(i).Transformation);
        
        % check for the personal best difference
        if particle(i).Distance < particle(i).Best.Distance
            particle(i).Best.Distance = particle(i).Distance;
            particle(i).Best.Transformation = particle(i).Transformation;
        end
        % check for new global best difference
        if particle(i).Best.Distance < GlobalBest.Distance
            GlobalBest = particle(i).Best;
        end
    end
    BestDistances(it) = GlobalBest.Distance;
    w = w*wdamp;
    disp(['Iteration ' num2str(it) ': Best Distance = ' num2str(BestDistances(it))]);
    
end
ObjectMoveNew = transformation(GlobalBest.Transformation, ObjectMove);
end