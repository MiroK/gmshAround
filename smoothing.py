from branching import find_branches, Walker, vertex_to_DG0_foo, smooth_data
from dolfin import MeshFunction, cells, Mesh, pi
import matplotlib.pyplot as plt
import numpy as np


mesh = Mesh('./data/vasc_mesh.xml.gz')
data = MeshFunction('int', mesh, './data/widths.xml.gz')

terminals, branches = find_branches(mesh)

# Smoothing
w = Walker(mesh)

plt.figure()

# Let's see how badly the data oscillates
for branch in branches:
    x, y = [], []
    for p, d in w.walk(branch):
        x.append(d)
        y.append(data[int(p)])
    plt.plot(x, y)

smooth_data(data)    

plt.figure()

# After smoothing
for branch in branches:
    x, y = [], []
    for p, d in w.walk(branch):
        x.append(d)
        y.append(data[int(p)])
    plt.plot(x, y)

# Conversion to something FEniCS_ii can use. NOTE that the data here
# should be on the 1d mesh that is found in 3d.
f = vertex_to_DG0_foo(data)
    
# Some statics
print 'Num branches', len(branches)

print 'Total length', sum(cell.volume() for cell in cells(mesh))
lengths = np.array([cell.volume() for cell in cells(mesh)])
branch_lengths = map(lambda branch: sum(lengths[c] for c in branch), branches)
print 'Longest/shortest branches', max(branch_lengths), min(branch_lengths)

branch_radii = f.vector().array()

dm = f.function_space().dofmap()
dm = [dm.cell_dofs(c.index())[0] for c in cells(mesh)]
branch_volumes = pi*branch_radii**2*lengths[dm]

print 'Total volume', sum(branch_volumes)
print 'Largest/smallest branches (by volume)', max(branch_volumes), min(branch_volumes)

plt.figure()
# Distribution of radii
n, bins, patches = plt.hist(branch_radii, facecolor='green')

plt.figure()
# Distribution of volume
n, bins, patches = plt.hist(branch_volumes, facecolor='blue')

plt.show()
# FIXME: is it a fractal
