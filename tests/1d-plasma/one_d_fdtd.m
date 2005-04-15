% 1d FDTD Simulation to try to determine why my 3d plasma is
% unstable.

len = 149;

c = 3e8;

centre_f = 500e12;
df = 200e12;

wl = c ./ (centre_f + df);

dz = wl ./ 32;
dt = 0.9 .* 1 ./ (c * sqrt(1/(dz*dz)));

disp(sprintf('dz = %g, dt = %g', dz, dt));
%dz = 18.75e-9;
%dt = dz / (2 * c);
%dt = 3.57482745e-17;

ex = zeros(len,1);
hy = zeros(len,1);
dx = ex;
sx = ex;
sxm1 = ex;
sxm2 = ex;

material = ones(len, 1);
material(80:135) = 2;

eps = [1 1] * 8.854185336732027e-12;
mu = [1 1] * 1.256637061435917e-06;
sigma = [0 0];
sigmastar = [0 0];
%plasma_f = [0 1.85e+15];
%collision_f = [0 1.4e+14];
%plasma_f = [0 5.8e15];
%collision_f = [0 107.5e12];
%plasma_f = [0 12.57e15];
%collision_f = [0 385.1e12];

% Johnson and Christy silver
plasma_f = [0 2.17e15];
collision_f = [0 32.258e12];

%mat.set_collision_freq(5.7e13 * 2 * M_PI);
%mat.set_plasma_freq(2e15 * 2 * M_PI);

Ca = (1 - (sigma * dt * 0.5) ./ eps) ./ ...
     (1 + (sigma * dt * 0.5) ./ eps);

Da = (1 - (sigmastar * dt * 0.5) ./ mu) ./ ...
     (1 + (sigmastar * dt * 0.5) ./ mu);

Cbz = (dt ./ (eps .* dz)) ./ (1 + (sigma .* dt .* 0.5) ./ eps);
Dbz = (dt ./ (mu .* dz)) ./ (1 + (sigmastar .* dt .* 0.5) ./ mu);

A = exp(-1 .* collision_f .* dt);
B = [0 (plasma_f(2) .* 2 .* pi) .^ 2 .* (dt ./ collision_f(2))];

ex40 = [];
ex60 = [];

%DFTfreqs = [300e12:10e12:800e12];
DFTfreqs = [centre_f - df:centre_f:centre_f + df];
ex40DFT_r = zeros(length(DFTfreqs), 1);
ex40DFT_i = zeros(length(DFTfreqs), 1);
ex60DFT_r = zeros(length(DFTfreqs), 1);
ex60DFT_i = zeros(length(DFTfreqs), 1);

ex2 = [0 0 0];
ex98 = [0 0 0];

figure
for idx = [1:1000]
  hy = update_hy(len, material, hy, Da, Dbz, ex);
  
  [ex, dx, sx, sxm1, sxm2] = update_ex(len, material, ex, Ca, Cbz, ...
                                       hy, A, B, dx, sx, sxm1, sxm2);
  
  %temp = gaussm(idx, 500e12, 4000e12, dt);
  %temp = gaussm(idx, 200e12, 100e12, dt);
  temp = gaussm(idx, df, centre_f, dt);
  %temp = 100 * gauss(idx, 200e12, dt);
  ex(20) = ex(20) + temp;

  % Neumann boundary
  %ex(1) = ex(2);
  %ex(len) = ex(len-1);
  
  % Dirchlet boundary
  %ex(1) = 0;
  %ex(len) = 0;
  
  % Simple absorber
  ex2 = [ex(2) ex2(1) ex2(2)];
  ex98 = [ex(148) ex98(1) ex98(2)];
  ex(1) = ex2(3);
  ex(149) = ex98(3);
  
  plot(ex); %, [1:len], hy);
  pause(0.1);
  %legend('Ex', 'Hy');
  %pause(0.25);
  
  ex40 = [ex40 ex(40)];
  ex60 = [ex60 ex(140)];
  
  [ex40DFT_r, ex40DFT_i] = point_dft(ex40DFT_r, ex40DFT_i, ...
                                     ex(40), DFTfreqs, dt * idx);
  [ex60DFT_r, ex60DFT_i] = point_dft(ex60DFT_r, ex60DFT_i, ...
                                     ex(140), DFTfreqs, dt * idx);
  
  %plot(ex60)
end


plot(DFTfreqs, sqrt(ex40DFT_r .^ 2 + ex40DFT_i .^ 2));
