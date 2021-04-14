function [dahd] = directed_averaged_hausdorff_distance(X,Y,step)
%Calculates the directed hausdorff distance for X to Y
%   If the directed Hausdorff distance is zero all points of X lie on a
%   point in Y, the converse might not be true

dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    fprintf('All points must have the same dimension')
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
dahd = sum(dXY_all)/(dimX(1)/step);
end

