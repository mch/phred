% Test sin modulated gaussian functions

deltat = 3.5748e-17;  % Time step size
N = 1500;         % Number of time points

f0 = 100e12;      % Centre frequency
df = 200e12;   % Frequency spacing

alpha = 1;      % Scaling factor

% Time data:
t = [0:deltat:deltat*(N-1)];

% Modulated pulse
%z = exp(-1 .* ((t-4 ./ (pi*df)) * df*pi).^2) .* sin(2*pi.*t*f0);
z = alpha * exp(-1 * ((t - 4. / (pi .* df)) .* df .* pi).^2) .* ...
    sin(2. .* pi .* f0 .* (t - 4. / (pi .* df)));
za = exp(-1 * ((t - 4. / (pi .* df)) .* df .* pi).^2);
zb = sin(2. .* pi .* f0 .* (t - 4. / (pi .* df)));
z2 = gauspuls(t, f0, 1.7); % Bandwidth of about 1.7*f0, ~500 THz. 

figure;
plot(t, z, t, za, t, zb, t, z2);
title('Sin modulated Gaussian Function');
legend('Jans FDTD implementation', 'Gaussian', 'sine', 'MATLAB gauspuls function');

Z = fft(z) ./ N;
Z2 = fft(z2) ./ N;

df = 1/(N*deltat);
%f = [-(N/2*df):df:(N/2 - 1)*df];
f = [0:df:(N/2 - 1)*df];

figure;
plot(f, abs(Z(1:N/2)), f, abs(Z2(1:N/2)));
title('Magnitude of frequncy response');
xlabel('Frequency (Hz)');
ylabel('Power');
legend('Jans', 'gauspuls');
