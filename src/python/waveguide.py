from Phred import *

temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Gedney's UPML in an X band waveguide, on %s %s..." % (MPI_SIZE,  temp)

# Waveguide size and length
a = 0.229 # m
b = 0.102 # m

centre_f = 9e9 # Hz
wl = 2.997925e8 / centre_f

# Extent of region
xlen = 200
ylen = 20
zlen = 45

# Prefix for output files
output_prefix = "wg2_"

num_time_steps = 3000

fdtd = FDTD()
fdtd.set_grid_size(xlen, ylen, zlen)
fdtd.set_grid_deltas(0.0051, 0.0051, 0.0051)

# Add boundary conditions
front = Pml()
back = Pml()
others = Ewall()

front.set_thickness(4)
back.set_thickness(4)

fdtd.set_boundary(FRONT, front)
#fdtd.set_boundary(BACK, back)
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

fdtd.add_geometry(interior)

# Excitation
gm = Gaussm()
gm.set_parameters(1, 0.2e9, 9e9)

ex = Excitation(gm)
ex.set_soft(0)
ex.set_region(10, 10, ylen/2, ylen/2, zlen/2, zlen/2);
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
pdft.freq_start = 8e9
pdft.freq_stop = 12.5e9
pdft.num_freqs = 120
pdft.set_point(p)
fdtd.add_result("centre_dft", pdft)
fdtd.map_result_to_datawriter("centre_dft", "mdw")

centrep = PointResult()
centrep.set_point(p)
fdtd.add_result("centrept", centrep)
fdtd.map_result_to_datawriter("centrept", "mdw")

p2 = point();
p2.x = 190
p2.y = ylen / 2
p2.z = zlen / 2

frontdft = PointDFTResult()
frontdft.freq_start = 8e9
frontdft.freq_stop = 12.5e9
frontdft.num_freqs = 120
frontdft.set_point(p2)
fdtd.add_result("front_dft", frontdft)
fdtd.map_result_to_datawriter("front_dft", "mdw")

frontpt = PointResult()
frontpt.set_point(p2)
fdtd.add_result("frontpt", frontpt)
fdtd.map_result_to_datawriter("frontpt", "mdw")

# Source output
adw = AsciiDataWriter(MPI_RANK, MPI_SIZE)
adw.set_filename(output_prefix + "source.txt")
fdtd.add_datawriter("adw", adw)

srcoutput = SourceTimeResult(gm)
fdtd.add_result("src", srcoutput)
fdtd.map_result_to_datawriter("src", "adw")
fdtd.map_result_to_datawriter("src", "mdw")

srcdft = SourceDFTResult(gm, 8e9, 12.5e9, 120)
fdtd.add_result("srcdft", srcdft)
fdtd.map_result_to_datawriter("srcdft", "mdw")

# Execute
fdtd.set_time_steps(num_time_steps)
fdtd.run(MPI_RANK, MPI_SIZE)
