from Phred import *

temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Gedney's UPML in an X band waveguide, on %s %s..." % (MPI_SIZE,  temp)

# Waveguide size and length
a = 0.229 # m
b = 0.102 # m

centre_f = 15e9 # Hz
wl = 2.997925e8 / centre_f

# Extent of region
xlen = 500
ylen = 317
zlen = 79

# Prefix for output files
output_prefix = "optm_" + str(MPI_SIZE) + "_"

num_time_steps = 1500

fdtd = FDTD()
fdtd.set_grid_size(xlen, ylen, zlen)
fdtd.set_grid_deltas(0.0001, 0.0001, 0.0001)
print "FDTD grid initialized. Cell size is 0.1mm x 0.1mm x 0.1mm\n"
print "Region is %ix%ix%i. Time step size is %g.\n" % (xlen, ylen,
                                                       zlen, fdtd.get_time_delta())


# Add boundary conditions
front = Pml()
back = Pml()
others = Ewall()

front.set_thickness(4)
back.set_thickness(4)

fdtd.set_boundary(FRONT, front)
fdtd.set_boundary(BACK, back)
#fdtd.set_boundary(BACK, others)
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

wall = Box()
wall.set_region(0, xlen, ylen/2, ylen/2, 0, zlen)
wall.material_id = 0

fdtd.add_geometry(interior)
fdtd.add_geometry(wall)

# Aperature's
num_aperatures = 3
ap_sizes = [70, 60] # Centre ap, Those to either side of centre

ap_pos = [xlen / 2, xlen / 2 + 70, xlen / 2 - 70]

# Aperatures are defined by thier position along the x-axis and thier
# size They are all square, a cannot be smaller than one cell, and
# larger than ylen. The aperture sizes and positions are symettric
# about the middle one.

if (1 == num_aperatures%2): # odd
    pass
else:
    pass

aps = [Box()]
aps[0].set_region(ap_pos[0] - ap_sizes[0]/2,
                  ap_pos[0] + ap_sizes[0]/2,
                  ylen/2, ylen/2,
                  zlen/2 - ap_sizes[0]/2,
                  zlen/2 + ap_sizes[0]/2)
aps[0].material_id = 1

for aperture in aps:
    fdtd.add_geometry(aperture)

                  

# Excitation
gm = Gaussm()
gm.set_parameters(1, 3e9, 15e9)

#ex = Excitation(gm)
#ex.set_soft(0)
#ex.set_region(10, 10, ylen/2, ylen/2, zlen/2, zlen/2);
#ex.set_polarization(0.0, 1.0, 0.0)

#exps = ExpSine(centre_f)
ex = WaveguideExcitation(gm)
ex.set_soft(1)
ex.set_region(15, 15, 0, ylen/2, 0, zlen)
ex.set_mode(0, 0, 1)
ex.set_polarization(0.0, 0.0, 1.0)

fdtd.add_excitation("modgauss", ex)

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
pr2.set_plane(p, TOP)
pr2.set_field(EY)
fdtd.add_result("ey_xyplane", pr2)
fdtd.map_result_to_datawriter("ey_xyplane", "ncdw")

pr3 = PlaneResult()
pr3.set_plane(p, TOP)
pr3.set_field(HX)
fdtd.add_result("hx_xyplane", pr3)
fdtd.map_result_to_datawriter("hx_xyplane", "ncdw")

pr4 = PlaneResult()
pr4.set_plane(p, TOP)
pr4.set_field(HZ)
fdtd.add_result("hz_xyplane", pr4)
fdtd.map_result_to_datawriter("hz_xyplane", "ncdw")

pdft = PointDFTResult()
pdft.freq_start = 12e9
pdft.freq_stop = 18e9
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
frontdft.freq_start = 12e9
frontdft.freq_stop = 18e9
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
fdtd.map_result_to_datawriter("src", "mdw")

srcdft = SourceDFTResult(gm, 12e9, 18e9, 120)
fdtd.add_result("srcdft", srcdft)
fdtd.map_result_to_datawriter("srcdft", "mdw")
fdtd.map_result_to_datawriter("srcdft", "adw")

# Execute
fdtd.set_time_steps(num_time_steps)
fdtd.run(MPI_RANK, MPI_SIZE)
