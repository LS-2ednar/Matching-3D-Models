function [hd] = hausdorff_distance(X, Y)
%calculates the hausdorff distance between two 
%point sets with the same dimension
dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    fprintf('All points must have the same dimension')
end

%calculate the directed distance from X to Y
dXY_all = zeros(1, dimX(1));
for i=1:dimX(1)
    shortestdist = norm(X(i,:) - Y(1,:));
    for j=2:dimY(1)
        dist = norm(X(i,:) - Y(j,:));
        if dist < shortestdist
            shortestdist = dist;
        end
    end
    dXY_all(i) = shortestdist;
end
dXY = max(dXY_all);
    
%calculate the directed averaged distance from Y to X    
dYX_all = zeros(1, dimY(1));
for j=1:dimY(1)
    shortestdist = norm(X(1,:) - Y(j,:));
    for i=2:dimX(1)
        dist = norm(X(i,:) - Y(j,:));
        if dist < shortestdist
            shortestdist = dist;
        end
    end
    dYX_all(j) = shortestdist;
end
dYX = max(dYX_all);

%calculate the modified hausdorff distance
hd = max(dXY, dYX);

end

