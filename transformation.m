function [X_new] = transformation(parameters, X)
%Transform a point cloud in 3D
alpha = parameters(1);
beta = parameters(2);
gamma = parameters(3);
xt = parameters(4);
yt = parameters(5);
zt = parameters(6);
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

T = [r11, r12, r13, r14;...
    r21, r22, r23, r24;...
    r31, r32, r33, r34;...
    r41, r42, r43, r44];

extension = ones(length(X),1);
X_ext = horzcat(X, extension);
X_new_ext = T*X_ext';
X_new = X_new_ext(1:3,:)';
end

