temp = "node"
if (MPI_SIZE > 1):
    temp = "nodes"
    
print "Testing the performance of Berengers PML on %s %s..." % (MPI_SIZE,  temp)

xlen = 100
ylen = 50
zlen = 50

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

# Fill interior
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

# Excitation

