function [f] = gaussm(ts, df, f0, dt)
% Gaussian modulated sin wave
t = ts .* dt;
f = exp(-((t-3 / (pi*df)) .* df * pi) .^ 2) ...
    .* sin(2 * pi * f0 .* (t - 4 / (pi * df)));
