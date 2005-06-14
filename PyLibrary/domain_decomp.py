# This is a test program for my domain decomposition algorithm as
# implemented in MPISubdomain.cc.

def decomp(grid_size, np):
    """Performs domain decomposition,

grid_size is a tuple containing the size of the grid (x, y, z)
np is the number of processors available
"""
    n = 1 # Divisions along x
    m = 1 # Divisions along y
    p = 1 # Divisions along z

    sdx = grid_size[0]
    sdy = grid_size[1]
    sdz = grid_size[2]

    divided = False

    for i in range(np):
        if (sdx >= sdy and sdx >= sdz and (n+1)*m*p <= np):
            n = n + 1
            sdx = grid_size[0] / n
            divided = True
            print "divided along x, n = %i, sdx = %f, sdy = %f, sdz = %f" \
                  % (n, sdx, sdy, sdz)

        elif (sdy >= sdx and sdy >= sdz and n*(m+1)*p <= np):
            m = m + 1
            sdy = grid_size[1] / m
            divided = True
            print "divided along y, m = %i, sdx = %f, sdy = %f, sdz = %f" \
                  % (m, sdx, sdy, sdz)
      
        elif (sdz >= sdx and sdz >= sdy and n*m*(p+1) <= np):
            p = p + 1
            sdz = grid_size[2] / p
            divided = True
            print "divided along z, p = %i, sdx = %f, sdy = %f, sdz = %f" \
                  % (p, sdx, sdy, sdz)
            
        if (divided == False):
            if ((n+1)*m*p <= np):
                n = n + 1
                sdx = grid_size[0] / n
                
            elif (n*(m+1)*p <= np):
                m = m + 1
                sdy = grid_size[1] / m

            elif (n*m*(p+1) <= np):
                p = p + 1;
                sdz = grid_size[2] / p

        divided = False

    print "Decomposition for a grid size of (%i, %i, %i) among %i processors:" \
          % (grid_size[0], grid_size[1], grid_size[2], np)
    print " -> %i, %i, %i" % (n, m, p)

    return (n, m, p)
