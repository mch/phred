function [dn dm dp] = dpoint(n, m, p, dx, dy, dz)

dn = (n.^2 .* dy .* dz - m.^2 .* dx .* dz + m.*dx.*dz - p.^2 .* dx .* dy + p .* ...
      dx .* dy) ./ (n.^2 .* m .* p);

dm = (m.^2 .* dx .* dz - n.^2 .* dy .* dz + n.*dy.*dz - p.^2 .* dx .* dy + p .* ...
      dx .* dy) ./ (n .* m.^2 .* p);

dp = (p.^2 .* dy .* dx - n.^2 .* dy .* dz + m.*dx.*dz - m.^2 .* dx .* dz + n .* ...
      dx .* dz) ./ (n .* m .* p.^2);
