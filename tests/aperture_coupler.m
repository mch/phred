function [N, V] = aperture_coupler()
% Your assignment, should you choose to accept it:
% Design an aperture couple and optimize the crap out of it. 

% Constants
A = 15.799e-3; % Ku band waveguide broad-wall size
B = 7.899e-3; % Ku band waveguide narrow-wall size

C = 3e8; % Speed of light, approximatly

EPS0 = 8.854185336732027e-12;
MU0 = 1.256637061435917e-06;
global C, EPS0, MU0;
C = 3e8; % Speed of light, approximatly

EPS0 = 8.854185336732027e-12;
MU0 = 1.256637061435917e-06;


MAX_APERTURES = 22;
global MAX_APERTURES;
MAX_APERTURES = 22;

% Parameters
f1 = 12e9; % Starting frequency
f2 = 18e9; % End frequency
f0 = 15e9; % Centre frequency
l0 = C / f0;

r = 0.1; % Input reflection coefficient
c = 0.7079; % Voltage coupling factor, 3dB
k1 = f1*2*pi * sqrt(MU0*EPS0);
k2 = f2*2*pi * sqrt(MU0*EPS0);
k0 = f0*2*pi * sqrt(MU0*EPS0);
kc = pi/A;
beta1 = sqrt(k1^2 - kc^2);
beta2 = sqrt(k2^2 - kc^2);
beta0 = sqrt(k0^2 - kc^2);

% Guided wave lengths
l_g1 = 2*pi/beta1;
l_g2 = 2*pi/beta2;
l_g0 = 2*pi/beta0;

% Coupler bandwidth
w_c = 2 * (l_g1 - l_g2) / (l_g1 + l_g2);

% Prototype transformer bandwidth
w_q = 2 * w_c;

% Max allowable VSWR
s_max = (1 + r) / (1 - r);

[N, V] = protyp(s_max, c, w_q);
Bn = sqrt(V);
Bn = Bn - 1./Bn % Normalized susceptances of our apertures

% Electrical lengths; distance from centre of one aperture to the
% next. 
phi_n = pi - 0.5 * (atan(2./Bn(1:length(Bn)-1)) + atan(2./Bn(2:length(Bn))))

% Optimization time! Find apertures which satisfy the Bn's. 
B_temp = Bn(1);
x_temp = 1:1:20;
y_temp = [];

for xt = x_temp
  y_temp(xt) = polarizability(B/xt, B, f0, l_g0, l0, 0.1e-3);
end

plot(x_temp, y_temp)




function [N, V] = protyp(s_max, c, w_q)
% A.3.6.2 converted to matlab
global MAX_APERTURES;
X0 = 1 / cos(0.25*pi*(2-w_q));
N = 1;
TN = zeros(MAX_APERTURES, 1);
TN(1) = X0;

psi = 2 * asin(c);

while (1)
  N = N + 1;
  
  if (N == MAX_APERTURES)
    break;
  end  

  psin = psi/N;
  S = (1 + sin(psin))/(1 - sin(psin));
  Sln = N * log(S);
  TN1 = TN(N-1);
  if (N == 2)
    TN2 = 1;
  else
    TN2 = TN(N-2);
  end
  TN(N) = 2 * X0 * TN1 - TN2;
  Ref = 1 + Sln / TN(N);
  if (Ref <= s_max)
    break;
  end
end

% Now, the rest of A.3.6.1 converted to matlab...
X = zeros(MAX_APERTURES);
FAK = zeros(MAX_APERTURES, 1);
V = zeros(MAX_APERTURES, 1);

X(1,1) = 2;
X(2,2) = X0;

if (N > 2)
  for I=3:N
    JU = 1;
    if (mod(I, 2) == 0)
      JU = 2;
    end
    
    for J = JU:2:I
      if (J == 1)
        X(I, J) = 2 * X0 * X(I-1, J+1) - X(I-2, J);
      else
        X(I, J) = X0 * (X(I-1, J-1) + X(I-1, J+1)) - X(I-2, J);
      end
    end
  end
end

N2 = (N+1)/2;

for I = 1:N2
  J = N + 2 - 2*I;
  FAK(I) = X(N,J);
  FAK(N+1 - I) = FAK(I);
end

FAC = FAK(1);
SUM = 0;

for I = 1:N
  FAK(I) = FAK(I) / FAC;
  SUM = SUM + FAK(I);
end

EXPO = Sln/SUM;

for I = 1:N
  V(I) = exp(FAK(I)*EXPO);
end



function bn = polarizability(a, b, f, l_g, l_0, t) % ? 
% Polarizability function for a rectangular aperture of width a and
% height b, from Dr. B's book.
% 
% Parameters: a, b, width and height of aperture, f is midband
% frequency, l_g is the guided wavelength, l_0 is the freespace
% wavelength, and t is the thickness of the wall in which the
% apertures are placed. 
%
% Returns the normalized susceptance produced by a given
% rectangular aperture
global C;

Fce = (C / 2) * sqrt((1/a^2) + (1/b^2));
Fcm = C / (2*a); 

Rm = 0;

if (b/a <= 0.5)
  Rm = (0.16695 * (b/a) + 0.0161)^0.8;
else
  Rm = 0.1939 * (b/a) + 0.061;
end

x = 0;

if (b <= a)
  x = b/a;
else
  x = a/b;
end

Re = 0.0925 * x^2 + 0.0209*x;

Am = 1.01095;
Ae = 0.45;

alpha_m = ((2*pi)/C) * sqrt(Fcm^2 - f^2);
alpha_e = ((2*pi)/C) * sqrt(Fce^2 - f^2);

Pm = (tan(pi*f/(2*Fcm)))/(pi*f/(2*Fcm)) * ...
     Rm * a^3 * exp(-alpha_m * Am * t)

Pe = (tan(pi*f/(2*Fce)))/(pi*f/(2*Fce)) * ...
     Re * a^3 * exp(-alpha_e * Ae * t)

% Position of the aperture along the x axis; x = 0 or A for an
% aperture in the narrow wall.
x = 0;

xx = (4*pi*Pm)/(a*b*l_g) * sin(pi*x*a)^2;

bz = (-1) * (pi*l_g*Pm)/(a^3*b) * cos(pi*x*a)^2;

by = (4*pi*l_g*Pe)/(a*b*l_0^2) * sin(pi*x*a)^2;

bn = 2 * abs(bz - xx + by);
