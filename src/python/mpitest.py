from Numeric import *

PI = 3.14159;

class ExpSource(SourceFunction):
    ampl_ = 1
    omega_ = 15e9 * 2 * PI;
    period_ = 1 / 15e9;
    def source_function(self, grid, time_step):
        t = grid.get_deltat() * time_step
        return self.ampl_ * (1 - exp(-1 * t / self.period_)) * sin(self.omega_ * t)


def mpitest(propagation_axis):
    temp = "node"
    if (MPI_SIZE > 1):
        temp = "nodes"
    
    print "Testing MPI in an X band waveguide, on %s %s..." % (MPI_SIZE,  temp)

    # Waveguide size and length
    a = 0.229 # m
    b = 0.102 # m

    centre_f = 15e9 # Hz
    wl = 2.997925e8 / centre_f

    # cells per wavelenth
    cpw = 15
    
    # cell size
    cs = wl / cpw;
    
    # Extent of region
    if (propagation_axis == 'x'):
        xlen = 270          
        ylen = int(a / cs)
        zlen = int(b / cs)
    elif (propagation_axis == 'y'):
        ylen = 270          
        xlen = int(a / cs)
        zlen = int(b / cs)
    elif (propagation_axis == 'z'):
        zlen = 270          
        ylen = int(a / cs)
        xlen = int(b / cs)
    else:
        print "Propagation axis: " + propagation_axis
        raise Exception("Invalid propagation axis given! Was given " )

    # Prefix for output files
    output_prefix = "mpi_test_" + str(MPI_SIZE) + "_" + str(xlen) + "_"

    num_time_steps = 110
    #num_time_steps = 500
    #num_time_steps = 10
    #num_time_steps = 270

    fdtd = FDTD()
    fdtd.set_grid_size(xlen, ylen, zlen)
    fdtd.set_grid_deltas(cs, cs, cs)
    
    print "FDTD grid initialized. Cell size is 0.1mm x 0.1mm x 0.1mm\n"
    print "Region is %ix%ix%i. Time step size is %g.\n" % (xlen, ylen,
                                                           zlen, fdtd.get_time_delta())
    print "Cells are %fx%fx%f m in size. " % (cs, cs, cs)

    # Add boundary conditions
    others = Ewall()
    
    fdtd.set_boundary(FRONT, others)
    fdtd.set_boundary(BACK, others)
    fdtd.set_boundary(LEFT, others)
    fdtd.set_boundary(RIGHT, others)
    fdtd.set_boundary(TOP, others)
    fdtd.set_boundary(BOTTOM, others)
    
    # Materials
    mlib = MaterialLib()
    mat = Material()
    mlib.add_material(mat)
    
    mat.epsilon = 2.2
    mat.name = "dielectric"
    mlib.add_material(mat)
    
    fdtd.load_materials(mlib)
    
    # Regions
    interior = Box()
    interior.set_region(0, xlen, 0, ylen, 0, zlen)
    interior.material_id = 1
    
    fdtd.add_geometry(interior)
    
    # Excitation
    #gm = Gaussm()
    #gm.set_parameters(1, 3e9, 15e9)
    
    #ex = Excitation(gm)
    #ex.set_soft(0)
    #ex.set_region(10, 10, ylen/2, ylen/2, zlen/2, zlen/2);
    #ex.set_polarization(0.0, 1.0, 0.0)

    #gm = ExpSine(centre_f)

    gm = ExpSource()
    
    #ex = WaveguideExcitation(gm)
    #ex.set_region(xlen/2, xlen/2, 0, ylen, 0, zlen)
    #ex.set_mode(0, 1, 0)

    ex = Excitation(gm)
    ex.set_region(xlen/2, xlen/2, ylen/4, ylen/4, zlen/2, zlen/2)

    ex.set_soft(1)
    ex.set_polarization(0.0, 0.0, 1.0)
    
    fdtd.add_excitation("modgauss", ex)
    
    # Data Writers!
    mdw = MatlabDataWriter(MPI_RANK, MPI_SIZE);
    mdw.set_filename(output_prefix + "point_data.mat")
    
    fdtd.add_datawriter("mdw", mdw)
    
    
    # Results!
    p = point();
    p.x = xlen / 2
    p.y = ylen / 2
    p.z = zlen / 2
    
    p1 = point();
    p1.x = (xlen / 2) - 50;
    p1.y = ylen / 2
    p1.z = zlen / 2
    
    p2 = point();
    p2.x = (xlen / 2) + 50
    p2.y = ylen / 2
    p2.z = zlen / 2
    
    print "Measurement point 1: %ix%ix%i, point 2: %ix%ix%i. " % (p1.x, p1.y, p1.z, p2.x, p2.y, p2.z)
    
    try:
        ncdw = NetCDFDataWriter(MPI_RANK, MPI_SIZE)
        ncdw.set_filename(output_prefix + "plane_data.nc")
        fdtd.add_datawriter("ncdw", ncdw)
        
        field_names = {EX: 'ex', EY: 'ey', EZ: 'ez', HX: 'hx', HY: 'hy', HZ: 'hz'}
        axis_faces = {'xy': TOP, 'yz': LEFT, 'xz': FRONT}
        for axis in ['xy', 'yz', 'xz']:
            for field in [EX, EY, EZ, HX, HY, HZ]:
                pr = PlaneResult()
                pr.set_name(field_names[field] + "_" + axis + "_plane")
                pr.set_field(field)
                pr.set_plane(p, axis_faces[axis])
                fdtd.add_result(field_names[field] + "_" + axis + "_plane", pr)
                fdtd.map_result_to_datawriter(field_names[field] + "_" + axis + "_plane", "ncdw")

    except Exception, e:
        print "Exception setting up NetCDF output:"
        print e

    p1r = PointResult()
    p1r.set_point(p1)
    fdtd.add_result("point1", p1r)
    fdtd.map_result_to_datawriter("point1", "mdw")

    p2r = PointResult()
    p2r.set_point(p2)
    fdtd.add_result("point2", p2r)
    fdtd.map_result_to_datawriter("point2", "mdw")
    
    # Source output
    srcoutput = SourceTimeResult(gm)
    fdtd.add_result("src", srcoutput)
    fdtd.map_result_to_datawriter("src", "mdw")
    
    srcdft = SourceDFTResult(gm, 12e9, 18e9, 120)
    fdtd.add_result("srcdft", srcdft)
    fdtd.map_result_to_datawriter("srcdft", "mdw")
    
    # Power output for S parameter measurement
    s11_r = region()
    s12_r = region()
    
    s11_r.xmin = p1.x; s11_r.xmax = p1.x;
    s11_r.ymin = 0; s11_r.ymax = ylen;
    s11_r.zmin = 0; s11_r.zmax = zlen;
    
    s12_r.xmin = p2.x; s12_r.xmax = p2.x;
    s12_r.ymin = 0; s12_r.ymax = ylen;
    s12_r.zmin = 0; s12_r.zmax = zlen;
    
    s11 = PowerResult(12e9, 18e9, 12)
    s11.set_region(s11_r)
    s12 = PowerResult(12e9, 18e9, 12)
    s12.set_region(s12_r)
    
    fdtd.add_result("s11", s11)
    fdtd.add_result("s12", s12)
    
    fdtd.map_result_to_datawriter("s11", "mdw")
    fdtd.map_result_to_datawriter("s12", "mdw")

    # Execute
    fdtd.set_time_steps(num_time_steps)
    fdtd.run(MPI_RANK, MPI_SIZE)


if (__name__ == "__main__"):
    mpitest('x')
   
