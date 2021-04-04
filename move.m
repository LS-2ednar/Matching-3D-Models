function [X] = move(X,Y)
%moves the first point set X to the center of gravity of the second point set

dimX = size(X);
dimY = size(Y);
%check if the two sets have the same dimension
if dimX(2) ~= dimY(2)
    fprintf('All points must have the same dimension')
end


%calculate the center of gravity for both point sets
center_of_gravity_X = zeros(1,dimX(2));
center_of_gravity_Y = zeros(1,dimY(2));

for i=1:dimX(2)
    x = sum(X(:,i))/dimX(1);
    y = sum(Y(:,i))/dimY(1);
    center_of_gravity_X(i) = x;
    center_of_gravity_Y(i) = y;
end

dist_vector = abs(center_of_gravity_X - center_of_gravity_Y);

for i=1:dimX(2)
    X(:,i) = X(:,i) + dist_vector(i);
end


end

