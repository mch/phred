function [C, Ceq] = nonlincon(X)
% Nonlinear constraints for the problem
C = [0; 0; 0; 0];
Ceq = [X(1)*X(2)*X(3) - 2];

% Force the X vector to contain only integers
Ceq = [Ceq; (X(1) - round(X(1)))^2; (X(2) - round(X(2)))^2; ...
       (X(3) - round(X(3)))^2];
