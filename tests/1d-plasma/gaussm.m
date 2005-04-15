function [f] = gaussm(ts, df, f0, dt)
% Gaussian modulated sin wave
t = ts .* dt;

%f = exp(-((t-3 / (pi*df)) .* df * pi) .^ 2) ...
%    .* sin(2 * pi * f0 .* (t - 4 / (pi * df)));

temp = (t - 4. / (pi * df)) * df * pi;

f = exp((-1) * temp * temp) ...
    * sin(2. * pi * f0 * (t - 4. / (pi * df)));
