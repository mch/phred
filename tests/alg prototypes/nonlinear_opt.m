% Do a little nonlinear optimization to see if this is gonna work
% or not. 

% Global domain dimensions:
dx = 100;
dy = 200;
dz = 1000;

% Number of processors
P = 2;

% Inital guess
n = 2;  % Number of subdomains along x axis
m = 1;  % Number of subdomains along y axis
p = 1;  % Number of subdomains along z axis

X = [n;m;p];

% Calculate point:
V = (n^2 * dy * dz - n * dy * dz + m^2 * dx * dz - m * dx * dz + p^2 ...
    * dx * dz - p * dx * dy) / (n*m*p);

% Gradient:

% Plot a few points
plot([1, 2, 3, 4, 5, 6, 7], [point([2; 1; 1]) ...
     point([1; 2; 1]) point([1; 1; 2]) point([2; 2; 1]) ...
     point([1; 2; 2]) point([2; 1; 2]) point([2;2;2])]);

points = [2 1 1; 1 2 1; 1 1 2; 2 2 1; 1 2 2; 2 1 2; 2 2 2];
% [dn, dm, dp] = dpoint(points(:,1), points(

% Nonlinear optimization
% The domain dimensions are set inside the point function
% definition, and the number of processors in set inside the
% nonlincon function definition. 

% Make n+m+p be greater than or equal to 2 (otherwise, whats the point?)
A = [1 1 1];
B = [4];
options = optimset('tolcon', 1e-7);
res = fmincon(@point, [n;m;p], A, B, [], [], [1;1;1], [],@nonlincon)
