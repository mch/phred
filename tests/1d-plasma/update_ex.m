function [ex, dx, sx, sxm1, sxm2] = update_ex(len, material, ex, Ca, ...
                                              Cbz, hy, A, B, dx, sx, ...
                                              sxm1, sxm2)

for j = [2:len]
  mid = material(j);
  
  if (A(mid) ~= 1)
    ex(j) = Ca(mid) * dx(j) + Cbz(mid) * (hy(j-1) - hy(j));
    dx(j) = ex(j);
    
    ex(j) = dx(j) - sx(j);

    sx(j) = (1 + A(mid)) * sxm1(j) - (A(mid) * sxm2(j)) + (B(mid) * ...
                                                  (1 - A(mid)) * ...
                                                  ex(j));
    sxm2(j) = sxm1(j);
    sxm1(j) = sx(j);
  else
    ex(j) = Ca(mid) * ex(j) + Cbz(mid) * (hy(j-1) - hy(j));
  end
end


