from Phred import *

temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Berengers PML on %s %s..." % (MPI_SIZE,  temp)

# Extent of region
xlen = 20
ylen = 15
zlen = 10

# Prefix for output files
output_prefix = "bpml_"

num_time_steps = 10

fdtd = FDTD()
fdtd.set_grid_size(xlen, ylen, zlen)
fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9)

# Add boundary conditions
left = Pml()
right = Pml()
top = Pml()
bottom = Pml()

ewall = Ewall()

left.set_thickness(4)
right.set_thickness(4)
top.set_thickness(4)
bottom.set_thickness(4)

fdtd.set_boundary(FRONT, ewall)
fdtd.set_boundary(BACK, ewall)
fdtd.set_boundary(LEFT, left)
fdtd.set_boundary(RIGHT, right)
fdtd.set_boundary(TOP, top)
fdtd.set_boundary(BOTTOM, bottom)

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
gm = Gaussm()
gm.set_parameters(10, 200e12, 100e12)

ex = Excitation(gm)
ex.set_soft(0)
ex.set_region(10, 10, ylen/2, ylen/2, zlen/2, zlen/2);
ex.set_polarization(0.0, 1.0, 0.0)

fdtd.add_excitation("modgauss", ex)

# Data Writers!
mdw = MatlabDataWriter(MPI_RANK, MPI_SIZE);
mdw.set_filename(output_prefix + "point_data_" + str(MPI_SIZE) + ".mat")

fdtd.add_datawriter("mdw", mdw)

# Results!
p = point();
p.x = xlen / 2
p.y = ylen / 2
p.z = zlen / 2

try:
    ncdw = NetCDFDataWriter(MPI_RANK, MPI_SIZE)
    ncdw.set_filename(output_prefix + "plane_data_" + str(MPI_SIZE) + ".nc")
    fdtd.add_datawriter("ncdw", ncdw)

    pr = PlaneResult()
    pr.set_plane(p, TOP)
    pr.set_field(EY)
    fdtd.add_result("xey", pr)
    fdtd.map_result_to_datawriter("xey", "ncdw")
    
    pr1 = PlaneResult()
    pr1.set_plane(p, FRONT)
    pr1.set_field(EX)
    fdtd.add_result("ex", pr1)
    fdtd.map_result_to_datawriter("ex", "ncdw")
    
    pr2 = PlaneResult()
    pr2.set_plane(p, FRONT)
    pr2.set_field(EY)
    fdtd.add_result("ey", pr2)
    fdtd.map_result_to_datawriter("ey", "ncdw")
    
    pr3 = PlaneResult()
    pr3.set_plane(p, FRONT)
    pr3.set_field(EZ)
    fdtd.add_result("ez", pr3)
    fdtd.map_result_to_datawriter("ez", "ncdw")
    
    pr1h = PlaneResult()
    pr1h.set_plane(p, FRONT)
    pr1h.set_field(HX)
    fdtd.add_result("hx", pr1h)
    fdtd.map_result_to_datawriter("hx", "ncdw")
    
    pr2h = PlaneResult()
    pr2h.set_plane(p, FRONT)
    pr2h.set_field(HY)
    fdtd.add_result("hy", pr2h)
    fdtd.map_result_to_datawriter("hy", "ncdw")
    
    pr3h = PlaneResult()
    pr3h.set_plane(p, FRONT)
    pr3h.set_field(HZ)
    fdtd.add_result("hz", pr3h)
    fdtd.map_result_to_datawriter("hz", "ncdw")
except:
    pass

pntr = PointResult()
pntr.set_point(p);
fdtd.add_result("pnt", pntr)
fdtd.map_result_to_datawriter("pnt", "mdw")

srcoutput = SourceTimeResult(gm)
fdtd.add_result("src", srcoutput)
fdtd.map_result_to_datawriter("src", "mdw")

print "p: (%i, %i, %i)" % (p.x, p.y, p.z)

# Execute
fdtd.set_time_steps(num_time_steps)
fdtd.run(MPI_RANK, MPI_SIZE)
