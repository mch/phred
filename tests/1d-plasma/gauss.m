function [f] = gauss(ts, df, dt)
% Gaussian modulated sin wave
t = ts .* dt;
f = exp(-((t-3 / (pi*df)) * df * pi) ^ 2);
