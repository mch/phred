function V = point(X)
% Calculates a point of the function we are minimizing. 
n = X(1);
m = X(2);
p = X(3);

% Region dimensions
dx = 100;
dy = 200;
dz = 1000;

%V = (n^2 * dy * dz - n * dy * dz + m^2 * dx * dz - m * dx * dz + p^2 ...
%    * dx * dz - p * dx * dy) / (n*m*p);

V = (n-1)*(dy*dz)/(m*p) + (m-1)*(dx*dz)/(n*p) + (p-1)*(dx*dy)/(n*m);
