print "Testing diffraction through a single hole in PEC...."

prefix = "hole_" + str(MPI_SIZE) + "_"
time_steps = 100
delta = 10e-9

fdtd = FDTD()
fdtd.set_grid_deltas(delta, delta, delta)
fdtd.set_grid_size(800e-9, 700e-9, 800e-9)

# Materials
mlib = MaterialLib()

mat = Material()
mat.epsilon = 2.2
mlib.add_material("dielectric", mat)

fdtd.load_materials(mlib)

# Boundaries
boundaries = [FRONT, BACK, LEFT, RIGHT, TOP, BOTTOM]

for b in boundaries:
    bound = Pml()
    bound.set_thickness(4)
    fdtd.set_boundary(b, bound)

ncdw = NetCDFDataWriter()
ncdw.set_filename(prefix + "planes.nc")
fdtd.add_datawriter("ncdw", ncdw)

mdw = MatlabDataWriter()
mdw.set_filename(prefix + "power.mat")
fdtd.add_datawriter("mdw", mdw)

#gridr = GridResult()
#fdtd.add_result("grid", gridr)
#fdtd.map_result_to_datawriter("grid", "ncdw")

fields = {FC_EZ: "ez", FC_EY: "ey", FC_EX: "ex", FC_HX: "hx", 
       FC_HY: "hy", FC_HZ: "hz"}

for f in fields.keys():
    pr = PlaneResult()
    pr.set_plane(grid_point(fdtd.get_x_cells() / 2, fdtd.get_y_cells() / 2, fdtd.get_z_cells() / 2), FRONT)
    pr.set_field(f)
    fdtd.add_result("plane_" + fields[f], pr)
    fdtd.map_result_to_datawriter("plane_" + fields[f], "ncdw")

#for f in fields.keys():
#    br = BlockResult()
#    br.set_time_param(10, 1000, 10)
#    br.set_field(f)
#    fdtd.add_result("block_" + fields[f], br)
#    fdtd.map_result_to_datawriter("block_" + fields[f], "ncdw")

splanes = CSGBox()
splanes.set_size(800e-9, 200e-9, 800e-9)

pinc = PowerResult(430e12, 600e13, 200)
pinc.set_region(splanes, LEFT)
fdtd.add_result("p_incident", pinc)

ptran = PowerResult(430e12, 600e13, 200)
ptran.set_region(splanes, RIGHT)
fdtd.add_result("p_transmission", ptran)

fdtd.map_result_to_datawriter("p_incident", "mdw")
fdtd.map_result_to_datawriter("p_transmission", "mdw")

gm = Gaussm()
gm.set_parameters(1, 86e12, 514e12)

st = SourceTimeResult(gm)
fdtd.add_result("src", st)
fdtd.map_result_to_datawriter("src", "mdw")

sdft = SourceDFTResult(gm, 430e12, 600e12, 200)
fdtd.add_result("srcdft", sdft)
fdtd.map_result_to_datawriter("srcdft", "mdw")

#ex = BartlettExcitation(gm)
ex = GaussWindow(gm)
exr = CSGBox()
exr.set_size(800e-9, delta * 2, 800e-9)
exr.set_centre(0, 250e-9, 0)
ex.set_region(exr)
ex.set_soft(1)
ex.set_type(E_FIELD)
ex.set_polarization(0, 0, 1)

fdtd.add_excitation("fluffy", ex)

fdtd.add_object("dielectric", exr)

# Make a hole in a metal plane:
metal = CSGBox()
metal.set_size(850e-9, 100e-9, 850e-9)
hole = CSGCylinder()
hole.set_radius(100e-9)
hole.set_height(200e-9)
rhole = CSGTransform(hole)
rhole.set_rotation(point(0,0,0), point(1, 0, 0), 90)

plate = CSGDifference(metal, rhole)

fdtd.add_object("PEC", plate)

fdtd.set_time_steps(time_steps)
fdtd.run()
