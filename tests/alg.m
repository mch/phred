% A little algorithm to divide a domain among identical processors

sdx = 100;
dx = 100;
n = 1;
sdy = 600;
dy = 200;
m = 1;
sdz = 1000;
dz = 1000;
p = 1;

% Number of processors
P = 8;
i = 0;
d = 0;
while (i < P)
  if (sdx >= sdy && sdx >= sdz && (n+1)*m*p <= P)
    n = n + 1;
    sdx = dx / n;
    d = 1;
    disp(sprintf(['Subdivided x, sdx = %d, sdy = %d, sdz =' ...
                  ' %d, (n,m,p) = (%d, %d, %d)'], sdx, sdy, sdz, n, ...
                 m, p));
  
  elseif (sdy >= sdx && sdy >= sdz && n*(m+1)*p <= P)
    m = m + 1;
    sdy = dy / m;
    d = 1;
    disp(sprintf(['Subdivided y, sdx = %d, sdy = %d, sdz =' ...
                  ' %d, (n,m,p) = (%d, %d, %d)'], sdx, sdy, sdz, n, ...
                 m, p));    

  elseif (sdz >= sdx && sdz >= sdy && (p+1)*n*m <= P)
    p = p + 1;
    sdz = dz / p;
    d = 1;
    disp(sprintf(['Subdivided z, sdx = %d, sdy = %d, sdz =' ...
                  ' %d, (n,m,p) = (%d, %d, %d)'], sdx, sdy, sdz, n, ...
                 m, p));
  end
  
  % If a division hasn't been sucessful, just do one based on the
  % number of processors. 
  if (d == 0)
    if ((n+1)*m*p <= P)
      n = n + 1;
      sdx = dx / n;
    elseif (n*(m+1)*p <= P)
      m = m + 1;
      sdy = dy / m;
    elseif (n*m*(p+1) <= P)
      p = p + 1;
      sdz = dz / p;
    end
  end
  
  i = i + 1;
  d = 0;
end

if (n*m*p ~= P)
  disp(['WARNING: domain decomp did not result in one subdomain for' ...
        ' each processor!']);
end


disp(sprintf(['A domain that is %dx%d%d is divided among %d processors' ...
              ' thusly:'], dx, dy, dz));
disp(sprintf(['The x axis is divided into %d pieces, each of' ...
              ' %d cells.'], n, sdx));
disp(sprintf(['The y axis is divided into %d pieces, each of' ...
              ' %d cells.'], m, sdy));
disp(sprintf(['The z axis is divided into %d pieces, each of' ...
              ' %d cells.'], p, sdz));
              

