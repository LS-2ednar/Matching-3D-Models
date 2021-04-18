clear;clc;close all;
rng(128);
%% ---------------------------- ATTENTION --------------------------------
%{
The total runtime of the script for both algorithms is about 125min 
(for 5 restarts). The average runtime of the simulated annealing algorithm 
is 24min (for 5 restarts of quick and fine alignment), 
the average runtime for the particle swarm  optimization is 100min 
(for 5 restarts of quick and fine alignment). 
PLEASE ADJUST THE NUMBER OF RESTARTS BELOW FOR A QUICKER RUNTIME.
(Lowering the number of iterations inside the optimization algorithms or
the particle swarm size is also possible but will yield considerable worse
results.)
%}
% number of restarts
numberOfrestarts = 5;
%% --------------------------- INTRODUCTION ------------------------------
%{
A problem in facial surgery is the reconstruction of the lower jaw. Which, 
for example, was destroyed by tumors or an accident. 
To reconstruct the lower jaw, a piece of bone is removed from the pelvic 
bone. In the following, the two algorithms Simulated Anneling and Particle 
Swarm will be used to determine an optimal point for cutting out the 
mandibular bone fragment. The similarity of the point set is based on the 
calculation of the directed averaged Hausdorff distance. A statement is 
made about the more suitable algorithm for this type of problem.
%}
%% -------- SETUP: Load Data & Create and Visualize Point Clouds ---------
%{
For an initial overview of the data, we plot the two 3D point sets
To get a better starting point for the initialization the mandible point 
sets center of gravity is moved to the center of gravity of the pelvis 
point set
%} 
tic
% reading in the stl files
stlData = stlread('Mand-left-cut.stl');
mand = stlData.Points;
stlData1 = stlread('Pelvis-left-cut.stl');
pelvis = stlData1.Points;
% Initial overview 3D plots of both stl objects 
% figure(1),
% plot3(mand(:,1),mand(:,2),mand(:,3),'.');
% title('Initial Point Cloud: Mandible')
% figure(2),
% plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
% title('Inital Point Cloud: Pelvis')

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
disp(['ELAPSED TIME FOR SETUP: ' num2str(toc) 's'])

%% ------------------ OPTIMIZATION: SIMULATED ANNEALING ------------------
%{ 
The first step is a quick alignment, the temperature is set to 50 and the
directed averaged hausdorff distance is calculated based on every 10th
point in both point sets. The quick alignment is followed by a fine
alignment, where we take the position of the mandible calculated in the
step before and only a small subsection of the pelvis where the mandible
is located now and run the algorithm again. The maximum step and rotation
is chosen significantly smaller than before, allowing the mandible to
move only slightly per cooling step. Since the calculation of the directed 
averaged hausdorff distance is still very computational expensive we
based the calculation on every second point in both point sets to achieve
a better runtime. 
%} 

tic
% Quick Alignment Parameters
startT_Q = 50;               % start temperature for the outside loop 
nInnerLoop_Q = 5;            % iterations for inner loop
maxStep_Q = 5;               % maximum translation step
maxRotation_Q = pi;          % maximum rotation 
acceptance_constant_Q = 80;  % determines the probability of accepting a
                             % suboptimal solution - for the quick
                             % alignment we want a higher probability, 
                             % to escape local minima 
dahd_step_Q = 10;            % stepsize for the hausdorff distance

% Fine Alignment Parameters
startT_F = 20;               % start temperature for the outside loop 
nInnerLoop_F = 5;            % iterations for inner loop
maxStep_F = 0.1;             % maximum translation step
maxRotation_F = 0.1*pi;      % maximum rotation
acceptance_constant_F = 200; % determines the probability of accepting a
                             % suboptimal solution - for the fine tuning
                             % we want a smaller probability!
dahd_step_F = 2;             % stepsize for the hausdorff distance

% initialize the mandible structure to save the results for each restart
mandible.Points = mand;
mandible.Distance = inf;

% save the best mand and its distance of the quick alignment for each restart 
MandQuick_SA = repmat(mandible, numberOfrestarts,1); 

% save the best mand and its distance of the fine alignment for each restart 
MandFine_SA = repmat(mandible, numberOfrestarts,1); 

for i = 1:numberOfrestarts
    disp(['SA: restart ' num2str(i)])
    [mand_quick_SA, distances_quick_SA] = SimulatedAnnealingOpti(mand,pelvis,...
        startT_Q, nInnerLoop_Q, maxStep_Q,maxRotation_Q,...
        acceptance_constant_Q,dahd_step_Q);

    MandQuick_SA(i).Points = mand_quick_SA;
    MandQuick_SA(i).Distance = min(distances_quick_SA);
    
    figure(4),
    plot(distances_quick_SA, 'Linewidth', 1);
    xlabel('Temperature')
    ylabel('Best Distance')
    title('SA: Best Distance per Temperature for Quick Alignment')
    set(gca,'Xdir','reverse')
    hold on
    drawnow
    
    figure(5),
    hold on
    plot3(mand_quick_SA(:,1),mand_quick_SA(:,2),mand_quick_SA(:,3),'.')
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
    title('SA: Best Transformations for Quick Alignment')
    drawnow
    
    % fine alignment 
    % subset the pelvis point set, to a small environment of the quick
    % solution for a quicker calculation of the hausdorff distance
    
    % since the mandible is not fully aligned with the pelvis point set,
    % we need to give the mandible more room to move in the subsection of
    % the pelvis
    correction_constant = 10; 
    xmin = min(mand_quick_SA(:,1))-correction_constant;
    xmax = max(mand_quick_SA(:,1))+correction_constant;
    ymin = min(mand_quick_SA(:,2))-correction_constant;
    ymax = max(mand_quick_SA(:,2))+correction_constant;
    zmin = min(mand_quick_SA(:,3))-correction_constant;
    zmax = max(mand_quick_SA(:,3))+correction_constant;
    
    pelvis_small = pelvis(pelvis(:,1) > xmin & pelvis(:,1) < xmax & ...
        pelvis(:,2) > ymin & pelvis(:,2) < ymax & ...
        pelvis(:,3) > zmin & pelvis(:,3) < zmax, :);
    
    % run the simulated annealing algorithm for the new mandible point set
    % whe obtained by the quick alignment and on the subsection of the pelvis
    [mand_fine_SA, distances_fine_SA] = SimulatedAnnealingOpti(mand_quick_SA,...
        pelvis_small,startT_F,nInnerLoop_F,maxStep_F,maxRotation_F,...
        acceptance_constant_F,dahd_step_F);

    MandFine_SA(i).Points = mand_fine_SA;
    MandFine_SA(i).Distance = min(distances_fine_SA);
    
    figure(6),
    plot(distances_fine_SA, 'Linewidth', 1);
    xlabel('Temperature')
    ylabel('Best Distance')
    title('SA: Best Distance per Temperature for Fine Alignment')
    set(gca,'Xdir','reverse')
    hold on
    drawnow
    
    figure(7),
    hold on
    plot3(mand_fine_SA(:,1),mand_fine_SA(:,2),mand_fine_SA(:,3),'.')
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
    title('SA: Best Transformations for Fine Alignment')
    drawnow
    
end
hold off

disp(['ELAPSED TIME FOR SIMULATED ANNEALING: ' num2str(toc) 's (for ' num2str(numberOfrestarts) ' restarts)'])

%% -------------------- OPTIMIZATION: PARTICLE SWARM ---------------------
%{
We again chose to do the alignment in two steps, to make the two algorithms
comparable we chose the same steps size for the calculation directed averaged 
hausdorff distance. The calculation is based on on every 10th point for the
quick alignment and and on every 2nd point for the fine alignment. 
The quick alignment has in total 50 iterations and a swarm size of 20 
particles. The fine alignment has only 20 iterations and a swarm size of 10 
particles, to safe some computation time.
%}
tic
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
iter_F = 20;             % number of iterations
nPop_F = 10;             % Population Size (Swarm Size)
wdamp_F = 0.9;           % damping coefficient
c1_F = 1.1;              % personal acceleration coefficent
c2_F = 1.1;              % social acceleration coefficient
stepsize_F = 2;          % stepsize for the hausdorff distance

mandible.Points = mand;
mandible.Distance = inf;

% save the best mand and its distance of the quick alignment for each restart 
MandQuick = repmat(mandible, numberOfrestarts,1); 
% save the best mand and its distance of the fine alignment for each restart 
MandFine = repmat(mandible, numberOfrestarts,1);

for i=1:numberOfrestarts
    disp(['PSO: restart ' num2str(i)])
    % quick alignment
    [mand_quick, distances_quick] = ParticleSwarmOpti(MinRot_Q, MaxRot_Q, MinTrans_Q, ...
        MaxTrans_Q, iter_Q, nPop_Q, wdamp_Q, c1_Q, c2_Q, mand, pelvis, stepsize_Q);
    
    MandQuick(i).Points = mand_quick;
    MandQuick(i).Distance = min(distances_quick);
    
    figure(8),
    plot(distances_quick, 'Linewidth', 1);
    xlabel('Iteration')
    ylabel('Best Distance')
    title('PSO: Best Distance per Iteration for Quick Alignment')
    hold on
    drawnow
    
    figure(9),
    hold on
    plot3(mand_quick(:,1),mand_quick(:,2),mand_quick(:,3),'.')
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
    title('PSO: Best Transformations for Quick Alignment')
    drawnow
   
    % fine alignment 
    % subset the pelvis point set, to a small environment of the quick
    % solution for a quicker calculation of the hausdorff distance
    
    % since the mandible is not fully aligned with the pelvis point set,
    % we need to give the mandible more room to move 
    correction_constant = 20;
    xmin = min(mand_quick(:,1))-correction_constant;
    xmax = max(mand_quick(:,1))+correction_constant;
    ymin = min(mand_quick(:,2))-correction_constant;
    ymax = max(mand_quick(:,2))+correction_constant;
    zmin = min(mand_quick(:,3))-correction_constant;
    zmax = max(mand_quick(:,3))+correction_constant;
    
    pelvis_small = pelvis(pelvis(:,1) > xmin & pelvis(:,1) < xmax & ...
        pelvis(:,2) > ymin & pelvis(:,2) < ymax & ...
        pelvis(:,3) > zmin & pelvis(:,3) < zmax, :);
    
    % run the particle swarm algorithm for the new mandible point set,
    % whe obtained by the quick alignment, and on the subsection of the pelvis
    [mand_fine, distances_fine] = ParticleSwarmOpti(MinRot_F, MaxRot_F, MinTrans_F, ...
        MaxTrans_F, iter_F, nPop_F, wdamp_F, c1_F, c2_F, mand_quick, pelvis_small, stepsize_F);
    
    MandFine(i).Points = mand_fine;
    MandFine(i).Distance = min(distances_fine);
    
    figure(10),
    plot(distances_fine, 'Linewidth', 1);
    xlabel('Iteration')
    ylabel('Best Distance')
    ylim([0 5])
    title('PSO: Best Distance per Iteration for Fine Alignment')
    hold on
    drawnow
    
    figure(11),
    hold on
    plot3(mand_fine(:,1),mand_fine(:,2),mand_fine(:,3),'.')
    plot3(pelvis(:,1),pelvis(:,2),pelvis(:,3),'.');
    title('PSO: Best Transformations for Fine Alignment')
    drawnow

end
hold off
disp(['ELAPSED TIME FOR PARTICLE SWARM OPTIMIZATION: ' num2str(toc) 's (for ' num2str(numberOfrestarts) ' restarts)'])

%% ----------------------------- CONCLUSION ------------------------------
%{
The simulated annealing algorithm runs considerably faster and produces 
point sets with distances that are very close to the one produced by PSO. 
Since it accepts suboptimal results with a higher probability
when the temperature is still high it can escape local optima and might find
the global one. In this example the algorithm has an optimum for the
mandible point set if we align it with the spine. Since this is not very
practical for the real world application we would have to adjust the
input point cloud for the pelvis and exclude the spine. 
While the running time of the particle swarm optimization does take a lot
longer it yields some interesting results. Every particle can 
'communicate' with all the other particles of the swarm, hence it is more 
likely to get stuck in a local optimum which can be intersting for the
quick alignment. 
One of the challenges for both algorithms lies in the definition and
calculation of the similarity function. The averaged directed hasudorff
distance seems to be not the best choice, since it computational expensive
when calculated for every point and does not yield good results when 
calculated for only a subsection of the points.
This was also an obstacle for parameter tuning since we wanted to get
reasonable results while still keeping the runtime as short as possible. 
%}
%% ----------------------------- FUNCTIONS -------------------------------

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
%             translation parameters
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

% extend the point cloud with 
extension = ones(length(X),1);
X_ext = horzcat(X, extension);

% transform the point cloud
X_new_ext = T*X_ext';
X_new = X_new_ext(1:3,:)';
end

function [dahd] = directed_averaged_hausdorff_distance(X,Y, step)
% Calculates the directed hausdorff distance from X to Y 
% If the directed averaged Hausdorff distance is zero all points of X lie on a
% point in Y, this can only happen if X <= Y.

% INPUT:
% X:    a 3D point set
% Y:    a 3D point set 
% OUTPUT:
% dahd: the directed averaged hausdorff distance from X to Y

dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    disp('All points must have the same dimension')
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

function [ObjectMoveNew, BestDistances] = SimulatedAnnealingOpti(ObjectMove,ObjectFixed,...
    startT,nInnerLoop,maxStep,maxRotation, acceptance_constant, dahd_step)
% Function for a the simulated annealing optimization for registration of 
% two 3D point clouds - find the best transformation of ObjectMove which
% minimized the directed average hausdorff distance to ObjectFixed

% INPUT:
% ObjectMove:           object which is moved though the soultion space
% ObjectFixed:          object which we want ObjectMove to align to
% startT:               inital temperatur
% nInnerLoop:           number of iterations through the inner loop
% maxStep:              maximal step size
% maxRotation:          maximal rotation 
% acceptance_constant:  determines the probability to accept suboptimal solutions 
% dahd_step:            directed avarage housdorff skiped units

% OUTPUT:
% ObjectMoveNew:        the moved 3D point cloud
% BestDistances:        the best distance for each Temperature


% INITIALIZATION
% Initialize paramters Rotation matrix to unit matrix and translation vector 
% to zero vector
% inital transformation parameters
alpha = 0;
beta = 0;
gamma = 0;
xt = 0;
yt = 0;
zt = 0;
parameters_best = [alpha, beta, gamma, xt, yt, zt];
parameters_current = parameters_best;

% Use hausdorff distance of the inital positions as initial best value
distance_best = directed_averaged_hausdorff_distance(ObjectMove, ObjectFixed, dahd_step);

% Calculate boundaries for the solution space
% the solution space is a cuboid around the pelvis point set, where each
% point of the point set has at least a distance of 5 points to the
% boarder of the cuboid. This allows the mandible to have at least an
% offset of 5 points from each point of the pelvis point cloud.

offset_constant = 5;
x_max = max(ObjectFixed(:,1))+offset_constant;
x_min = min(ObjectFixed(:,1))-offset_constant;
y_max = max(ObjectFixed(:,2))+offset_constant;
y_min = min(ObjectFixed(:,2))-offset_constant;
z_max = max(ObjectFixed(:,3))+offset_constant;
z_min = min(ObjectFixed(:,3))-offset_constant;

BestDistances = zeros(startT, 1);

% OPTIMIZATION
fprintf(['Max T:',num2str(startT),'\n'])
    for T=startT:-1:1
        for v=1:nInnerLoop
            % randomly update parameters for rotation
            % dynamic stepsize: the change of parameter becomes potentially
            % smaller when the temperature cools down
            parameters_current(1) = parameters_best(1) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(2) = parameters_best(2) + (rand-0.5)*2*maxRotation*T/startT;
            parameters_current(3) = parameters_best(3) + (rand-0.5)*2*maxRotation*T/startT;

            % randomly update parameters for translation
            parameters_current(4) = parameters_best(4) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(5) = parameters_best(5) + (rand-0.5)*2*maxStep*T/startT;
            parameters_current(6) = parameters_best(6) + (rand-0.5)*2*maxStep*T/startT;

            % transform the mandibe point set
            mand_current = transformation(parameters_current, ObjectMove);

            % update the parameters as long as we are not in the solution space
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
                mand_current = transformation(parameters_current, ObjectMove);
            end
            % calculated the (modified) hausdorff distance for the transformed
            % mand matrix 
            distance_current = directed_averaged_hausdorff_distance(mand_current, ObjectFixed,dahd_step);
            difference = distance_current - distance_best;

            % if the new distance is smaller than the last distance accept the
            % solution
            if difference < 0
                parameters_best = parameters_current;
                distance_best = distance_current;

            % else if the new distance is not smaller than the last distance
            % accept the solution with a random probability 
            elseif (exp((-difference*acceptance_constant)/T) > rand)
                parameters_best = parameters_current;
                distance_best = distance_current;
            end
        end
        BestDistances(T) = distance_best;
        disp(['Temperature ' num2str(T) ': Best Distance = ' num2str(BestDistances(T))]);
    end
    ObjectMoveNew = transformation(parameters_best, ObjectMove);
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
% BestDistance:         the best distance for each iteration

% The Objective Function we want to minimize
% x: vector of transformation parameters
ObjFunc = @(x) directed_averaged_hausdorff_distance(transformation(x, ObjectMove), ObjectFixed, stepsize);

nParameters = 6;    % number of transformation parameters
w = 1;              % Inertia coefficient

% INITIALIZATION
% initialize the transformation parameters to zero
% and Velocity and Distance as empty arrays
inital_particle.Transformation = zeros(1,nParameters);
inital_particle.Velocity = [];
inital_particle.Distance = [];
inital_particle.Best.Transformation = [];
inital_particle.Best.Distance = [];

particle = repmat(inital_particle, nPopulation,1); % initialize the whole population

GlobalBest.Distance = inf;                         % initialize the worst case

for i=1:nPopulation
    
    % set random values for the rotation parameters between MinRotation and
    % MaxRotation
    particle(i).Transformation(:,1:3) = unifrnd(MinRotation, MaxRotation, [1 3]);
    % set random values for the translation parameters between
    % MinTranslation and MaxTranslation
    particle(i).Transformation(:,4:6) = unifrnd(MinTranslation, MaxTranslation, [1 3]);
    
    % initialize the velocity to zero
    particle(i).Velocity = zeros([1 nParameters]);
    
    % evaluate the objective function to get the distance
    particle(i).Distance = ObjFunc(particle(i).Transformation);
    
    % set the 
    particle(i).Best.Transformation = particle(i).Transformation;
    particle(i).Best.Distance = particle(i).Distance;
    
    if particle(i).Best.Distance < GlobalBest.Distance
        GlobalBest = particle(i).Best;
    end

end

% hold best global distance value for each iteration
BestDistances = zeros(iter, 1);

% OPTIMIZATION
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
    
    % safe the best distance for each iteration
    BestDistances(it) = GlobalBest.Distance;
    
    % update the inertia coefficient
    w = w*wdamp;
    disp(['Iteration ' num2str(it) ': Best Distance = ' num2str(BestDistances(it))]);
    
end
% return the transformed object
ObjectMoveNew = transformation(GlobalBest.Transformation, ObjectMove);
end