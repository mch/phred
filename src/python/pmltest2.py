from Numeric import *

PI = 3.14159;

class SinSource(SourceFunction):
    ampl_ = 1
    omega_ = 15e9 * 2 * PI;
    period_ = 1 / 15e9;
    def source_function(self, grid, time_step):
        t = grid.get_deltat() * time_step
        return self.ampl_ * sin(self.omega_ * t)


def pml_test(pml_thickness, xlen):
    temp = "node"
    if (MPI_SIZE > 1):
        temp = "nodes"
    
    print "Testing the performance of Berengers's PML in an X band waveguide, on %s %s..." % (MPI_SIZE,  temp)

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
    #xlen = 270
    ylen = int(a / cs)
    zlen = int(b / cs)

    # Prefix for output files
    output_prefix = "pml_test_" + str(MPI_SIZE) + "_" + str(xlen) + "_" + str(pml_thickness) + "_"

    num_time_steps = 270
    #num_time_steps = 10

    fdtd = FDTD()
    fdtd.set_grid_size(xlen, ylen, zlen)
    fdtd.set_grid_deltas(cs, cs, cs)
    print "FDTD grid initialized. Cell size is 0.1mm x 0.1mm x 0.1mm\n"
    print "Region is %ix%ix%i. Time step size is %g.\n" % (xlen, ylen,
                                                           zlen, fdtd.get_time_delta())
    print "Cells are %fx%fx%f m in size. " % (cs, cs, cs)
    print "PML thickness is %i cells. " % pml_thickness

    # Add boundary conditions
    front = Pml()
    back = Pml()
    others = Ewall()
    
    front.set_thickness(pml_thickness)
    back.set_thickness(pml_thickness)
    
    fdtd.set_boundary(FRONT, front)
    fdtd.set_boundary(BACK, back)
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

    gm = SinSource()
    
    ex = WaveguideExcitation(gm)
    ex.set_soft(1)
    ex.set_region(xlen/2, xlen/2, 0, ylen, 0, zlen)
    ex.set_mode(0, 1, 0)
    ex.set_polarization(0.0, 0.0, 1.0)
    
    fdtd.add_excitation("modgauss", ex)
    
    # Data Writers!
    mdw = MatlabDataWriter();
    mdw.set_filename(output_prefix + "point_data.mat")
    
    fdtd.add_datawriter("mdw", mdw)
    
    
    # Results!
    p = point();
    p.x = xlen / 2
    p.y = ylen / 2
    p.z = zlen / 2
    
    p1 = point();
    p1.x = (xlen / 2) - 40;
    p1.y = ylen / 2
    p1.z = zlen / 2
    
    p2 = point();
    p2.x = (xlen / 2) - 40
    p2.y = ylen / 2 - ylen / 4
    p2.z = zlen / 2

    p3 = point();
    p3.x = (xlen / 2) - 5;
    p3.y = ylen / 2
    p3.z = zlen / 2

    p4 = point();
    p4.x = xlen / 2
    p4.y = ylen / 2
    p4.z = 5;

    p5 = point();
    p5.x = xlen / 2
    p5.y = ylen / 2
    p5.z = zlen / 2 + 10;
    
    print "Measurement point 1: %ix%ix%i, point 2: %ix%ix%i. " % (p1.x, p1.y, p1.z, p2.x, p2.y, p2.z)
    
    try:
        ncdw = NetCDFDataWriter()
        ncdw.set_filename(output_prefix + "plane_data.nc")
        fdtd.add_datawriter("ncdw", ncdw)
        
        pr1 = PlaneResult()
        # #pr1.set_name("ey-xzplane")
        pr1.set_plane(p4, TOP)
        pr1.set_field(EZ)
        fdtd.add_result("ez_xyplane_p4", pr1)
        fdtd.map_result_to_datawriter("ez_xyplane_p4", "ncdw")
        
        pr2 = PlaneResult()
        # #pr2.set_name("ey-xzplane")
        pr2.set_plane(p5, TOP)
        pr2.set_field(EZ)
        fdtd.add_result("ez_xyplane_p5", pr2)
        fdtd.map_result_to_datawriter("ez_xyplane_p5", "ncdw")
        
        pr3 = PlaneResult()
        # #pr2.set_name("ey-xzplane")
        pr3.set_plane(p, TOP)
        pr3.set_field(EZ)
        fdtd.add_result("ez_xyplane_c", pr3)
        fdtd.map_result_to_datawriter("ez_xyplane_c", "ncdw")
        
        pr4 = PlaneResult()
        pr4.set_plane(p1, FRONT)
        pr4.set_field(EZ)
        fdtd.add_result("ez_yzplane", pr4)
        fdtd.map_result_to_datawriter("ez_yzplane", "ncdw")    

        pr5 = PlaneResult()
        pr5.set_plane(p1, LEFT)
        pr5.set_field(EZ)
        fdtd.add_result("ez_xzplane", pr5)
        fdtd.map_result_to_datawriter("ez_xzplane", "ncdw")    

    except Exception, e:
        print "Script caught exception: "
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
    pw = region()

    pw.xmin = p3.x; pw.xmax = p3.x;
    pw.ymin = 0; pw.ymax = ylen;
    pw.zmin = 0; pw.zmax = zlen;
    
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

    pwr = PowerResult(12e9, 18e9, 12)
    pwr.set_region(pw)
    
    fdtd.add_result("s11", s11)
    fdtd.add_result("s12", s12)
    fdtd.add_result("pwr", pwr)
    
    fdtd.map_result_to_datawriter("s11", "mdw")
    fdtd.map_result_to_datawriter("s12", "mdw")
    fdtd.map_result_to_datawriter("pwr", "mdw")

    # Execute
    fdtd.set_time_steps(num_time_steps)
    fdtd.run(MPI_RANK, MPI_SIZE)


if (__name__ == "__main__"):
    pml_test(4, 110)
    #pml_test(4, 165)
    #pml_test(8, 110)
    #pml_test(8, 165)
