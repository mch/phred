temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Berengers PML on %s %s..." % (MPI_SIZE,  temp)

# Extent of region
xlen = 100
ylen = 50
zlen = 50

# Centre of sphere:
centre = Types.point()
centre.x = 40
centre.y = ylen / 2
centre.z = zlen / 2

num_time_steps = 500

fdtd = FDTD.FDTD()
fdtd.set_grid_size(xlen, ylen, zlen)
fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9)

# Add boundary conditions
front = Boundaries.Pml()
back = Boundaries.Pml()
left = Boundaries.Pml()
right = Boundaries.Pml()
top = Boundaries.Pml()
bottom = Boundaries.Pml()

front.set_thickness(4)
back.set_thickness(4)
top.set_thickness(4)
bottom.set_thickness(4)
left.set_thickness(4)
right.set_thickness(4)

fdtd.set_boundary(Types.FRONT, front)
fdtd.set_boundary(Types.BACK, back)
fdtd.set_boundary(Types.LEFT, left)
fdtd.set_boundary(Types.RIGHT, right)
fdtd.set_boundary(Types.TOP, top)
fdtd.set_boundary(Types.BOTTOM, bottom)

# Materials
mlib = Materials.MaterialLib()
mat = Materials.Material()
mlib.add_material(mat)

mat.epsilon = 2.2
mat.name = "dielectric"
mlib.add_material(mat)

mat.epsilon = 1
mat.name = "silver"
mat.collision_freq = 1.4e14
mat.plasma_freq = 2 * 3.14159 * 1.85e15
mlib.add_material(mat)

fdtd.load_materials(mlib)

# Regions
interior = Geometry.Box()
interior.set_region(0, xlen, 0, ylen, 0, zlen)
interior.material_id = 1

sphere = Geometry.Sphere(centre, 5)
sphere.material_id = 3

fdtd.add_geometry(interior)
fdtd.add_geometry(sphere)

# Excitation
gm = Excitations.Gaussm()
gm.set_parameters(10, 200e12, 100e12)

ex = Excitations.Excitation(gm)
ex.set_soft(0)
ex.set_region(20, 20, ylen/2, ylen/2, zlen/2, zlen/2);
ex.set_polarization(0.0, 1.0, 0.0)

fdtd.add_e_excitation("modgauss", ex)

# Data Writers!
mdw = DataWriters.MatlabDataWriter(MPI_RANK, MPI_SIZE);
mdw.set_filename("test.mat")

fdtd.add_datawriter("mdw", mdw)

ncdw = DataWriters.NetCDFDataWriter(MPI_RANK, MPI_SIZE)
ncdw.set_filename("yz_plane.nc")
fdtd.add_datawriter("ncdw", ncdw)

# Results!
p = Types.point();
p.x = xlen / 2
p.y = ylen / 2
p.z = zlen / 2

pr2 = Results.PlaneResult()
#pr2.set_name("ey-xzplane")
pr2.set_plane(p, Types.LEFT)
pr2.set_field(Types.EY)

# Execute
fdtd.set_time_steps(num_time_steps)
fdtd.run(MPI_RANK, MPI_SIZE)
