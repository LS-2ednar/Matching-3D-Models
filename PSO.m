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
figure,
plot3(mand(:,1),mand(:,2),mand(:,3),'.')
hold on
plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'k.');
title('The inital position of the mand and pelvis')


%% Particle Swarm Optimization

CostFunction = @(x) directed_averaged_hausdorff_distance(transformation(x, mand), pelvis, 20);


nVar = 6; % number of parameters
VarSize = [1 nVar]; 

VarMinRotation = 0;
VarMaxRotation = 2*pi;

VarMinTranslation = -50;
VarMaxTranslation = 50;

%% Parameters

MaxIt = 100; % Number of Iteration
nPop = 20; % Population Size (Swarm Size)

w = 1; % Inertia coefficient
wdamp = 0.99; 
c1 = 1; % personal acceleration coefficent
c2 = 1; % social acceleration coefficient

%% Initialization

empty_particle.Position = zeros(1,6);
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];


particle = repmat(empty_particle, nPop,1);


GlobalBest.Cost = inf;


for i=1:nPop
    
    particle(i).Position(:,1:3) = unifrnd(VarMinRotation, VarMaxRotation, [1 3]);
    particle(i).Position(:,4:6) = unifrnd(VarMinTranslation, VarMaxTranslation, [1 3]);
    
    particle(i).Velocity = zeros(VarSize);
    
    particle(i).Cost = CostFunction(particle(i).Position);
    
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end

end

% Hold best cost value for each iteration
BestCosts = zeros(MaxIt, 1);


%% Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        particle(i).Velocity = w*particle(i).Velocity + ...
            c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) + ...
            c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        particle(i).Cost = CostFunction(particle(i).Position);
        
        if particle(i).Cost < particle(i).Best.Cost
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Position = particle(i).Position;
        end
        
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end
    end
    BestCosts(it) = GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    
    w = w*wdamp;
    
end

%% Results
mand = transformation(GlobalBest.Position, mand);
figure
mand_fig = plot3(mand(:,1),mand(:,2),mand(:,3),'.');
hold on
pelvis_fig = plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
title('Best Alignment')
