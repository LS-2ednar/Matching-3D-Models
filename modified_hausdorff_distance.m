function [mhd] = modified_hausdorff_distance(X,Y)
%calculates the modified hausdorff distance between two 
%point sets with the same dimension

dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    fprintf('All points must have the same dimension')
end

%calculate the directed averaged distance from X to Y
dXY = 0;
for i=1:dimX(1)
    shortestdist = norm(X(i,:) - Y(1,:));
    for j=2:dimY(1)
        dist = norm(X(i,:) - Y(j,:));
        if dist < shortestdist
            shortestdist = dist;
        end
    end
    dXY = dXY + shortestdist;
end
adXY = dXY/dimX(1);
    
%calculate the directed averaged distance from Y to X    
dYX = 0;
for j=1:dimY(1)
    shortestdist = norm(X(1,:) - Y(j,:));
    for i=2:dimX(1)
        dist = norm(X(i,:) - Y(j,:));
        if dist < shortestdist
            shortestdist = dist;
        end
    end
    dYX = dYX + shortestdist;
end
adYX = dYX/dimY(1);

%calculate the modified hausdorff distance
mhd = max(adXY, adYX);

end

