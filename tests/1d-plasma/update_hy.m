function hy = update_hy(len, material, hy, Da, Dbz, ex)

for j = [1:len-1]
  mid = material(j);
  hy(j) = Da(mid) * hy(j) + Dbz(mid) * (ex(j) - ex(j+1));
end
