from Phred import *

temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Gedney's UPML on %s %s..." % (MPI_SIZE,  temp)

# Extent of region
xlen = 100
ylen = 50
zlen = 50

# Centre of sphere:
centre = point()
centre.x = 40
centre.y = ylen / 2
centre.z = zlen / 2

# Prefix for output files
output_prefix = "upml_"

num_time_steps = 500

fdtd = FDTD()
fdtd.set_grid_size(xlen, ylen, zlen)
fdtd.set_grid_deltas(18.75e-9, 18.75e-9, 18.75e-9)

# Add boundary conditions
front = UPml()
back = UPml()
left = UPml()
right = UPml()
top = UPml()
bottom = UPml()

front.set_thickness(4)
back.set_thickness(4)
top.set_thickness(4)
bottom.set_thickness(4)
left.set_thickness(4)
right.set_thickness(4)

fdtd.set_boundary(FRONT, front)
fdtd.set_boundary(BACK, back)
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

mat.epsilon = 1
mat.name = "silver"
mat.collision_freq = 1.4e14
mat.plasma_freq = 2 * 3.14159 * 1.85e15
mlib.add_material(mat)

fdtd.load_materials(mlib)

# Regions
interior = Box()
interior.set_region(0, xlen, 0, ylen, 0, zlen)
interior.material_id = 1

sphere = Sphere(centre, 5)
sphere.material_id = 3

fdtd.add_geometry(interior)
#fdtd.add_geometry(sphere)

# Excitation
gm = Gaussm()
gm.set_parameters(10, 200e12, 100e12)

ex = Excitation(gm)
ex.set_soft(0)
ex.set_region(20, 20, ylen/2, ylen/2, zlen/2, zlen/2);
ex.set_polarization(0.0, 1.0, 0.0)

fdtd.add_e_excitation("modgauss", ex)

# Data Writers!
mdw = MatlabDataWriter(MPI_RANK, MPI_SIZE);
mdw.set_filename(output_prefix + "point_data.mat")

fdtd.add_datawriter("mdw", mdw)

ncdw = NetCDFDataWriter(MPI_RANK, MPI_SIZE)
ncdw.set_filename(output_prefix + "plane_data.nc")
fdtd.add_datawriter("ncdw", ncdw)

# Results!
p = point();
p.x = xlen / 2
p.y = ylen / 2
p.z = zlen / 2

pr2 = PlaneResult()
#pr2.set_name("ey-xzplane")
pr2.set_plane(p, LEFT)
pr2.set_field(EY)
fdtd.add_result("ey_xzplane", pr2)
fdtd.map_result_to_datawriter("ey_xzplane", "ncdw")

pdft = PointDFTResult()
pdft.freq_start = 5e12
pdft.freq_stop = 600e12
pdft.num_freqs = 120
pdft.set_point(p)
fdtd.add_result("centre_dft", pdft)
fdtd.map_result_to_datawriter("centre_dft", "mdw")

centrep = PointResult()
centrep.set_point(p)
fdtd.add_result("centrept", centrep)
fdtd.map_result_to_datawriter("centrept", "mdw")

p2 = point();
p2.x = 10
p2.y = ylen / 2
p2.z = zlen / 2

frontdft = PointDFTResult()
frontdft.freq_start = 5e12
frontdft.freq_stop = 600e12
frontdft.num_freqs = 120
frontdft.set_point(p2)
fdtd.add_result("front_dft", frontdft)
fdtd.map_result_to_datawriter("front_dft", "mdw")

frontpt = PointResult()
frontpt.set_point(p2)
fdtd.add_result("frontpt", frontpt)
fdtd.map_result_to_datawriter("frontpt", "mdw")


# Execute
fdtd.set_time_steps(num_time_steps)
fdtd.run(MPI_RANK, MPI_SIZE)
