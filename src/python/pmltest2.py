from Phred import *

temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Berengers PML on %s %s..." % (MPI_SIZE,  temp)

# Extent of region
xlen = 900
ylen = 50
zlen = 50

# Prefix for output files
output_prefix = "pmltest_" + str(xlen) + "_"

num_time_steps = 900

fdtd = FDTD()
fdtd.set_grid_size(xlen, ylen, zlen)
fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9)

# Add boundary conditions
top = Pml()

ewall = Ewall()

top.set_thickness(4)

fdtd.set_boundary(FRONT, ewall)
fdtd.set_boundary(BACK, ewall)
fdtd.set_boundary(LEFT, ewall)
fdtd.set_boundary(RIGHT, ewall)
fdtd.set_boundary(TOP, top)
fdtd.set_boundary(BOTTOM, ewall)

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
p.x = 90
p.y = ylen / 2
p.z = zlen / 2

pntr = PointResult()
pntr.set_point(p);
fdtd.add_result("pnt", pntr)
fdtd.map_result_to_datawriter("pnt", "mdw")

print "p: (%i, %i, %i)" % (p.x, p.y, p.z)

# Execute
fdtd.set_time_steps(num_time_steps)
fdtd.run(MPI_RANK, MPI_SIZE)
