# This script illustrates how the Xd-1d mesh could be makde from 1d
from dolfin import *
from random import sample
from fenics_ii.trace_tools.embedded_mesh import EmbeddedMesh
from mesh_around_1d import mesh_around_1d
from meshconvert import convert2xml
import subprocess, os
import numpy as np

# Fabricate some 1d mesh and store it together with some vertex function
mesh = UnitSquareMesh(32, 32)
mesh.init(1)

n = mesh.num_entities(1)
f = FacetFunction('size_t', mesh, 0)

for i in sample(range(n), n/4): f[i] = 1

# The 1d mesh
mesh = EmbeddedMesh(mesh, f, 1).mesh
# Fake data
data = VertexFunction('double', mesh, 0)

# Set to distance from ori so that we can check later
origin = Point(0., 0.)
for vertex in vertices(mesh):
    data[vertex] = vertex.point().distance(origin)

out, out_data = 'foo.xml', 'foo_data.xml'
File(out) << mesh
File(out_data) << data

# Create the 2d mesh around it
geo, _ = mesh_around_1d(out)

ccall= 'gmsh -%d -optimize %s' % (2, geo)
subprocess.call(ccall, shell=True)
    
# Convert
xml_file = 'foo.xml'
msh_file = 'foo.msh'
convert2xml(msh_file, xml_file)

# Throw away the line function as it is expensive to store all the edges
meshXd = Mesh(xml_file)
# Make sure that the assumption of correspondence hold
assert np.linalg.norm(mesh.coordinates() -
                      meshXd.coordinates()[:mesh.num_vertices()]) < 1E-13

line_f = MeshFunction('size_t', meshXd, 'foo_facet_region.xml')
assert all(line_f[e] == 1 or line_f[e] == 0 for e in edges(meshXd))

# Just a list of edges which are 1
np.savetxt('foo_1d.txt', [edge.index() for edge in SubsetIterator(line_f, 1)])


# Data from 1d can be represented as a list
np.savetxt('foo_1d_data.txt',
           [[vertex.index(), data[vertex]] for vertex in vertices(mesh)])
# Original 1d

################
# Loading part
################
# Now we going to build the 1d mesh back and get the data as well
mesh = Mesh('foo.xml')  # 2d

# FIXME: embedded mesh should accept a list
embedded_mark = EdgeFunction('size_t', mesh, 0)
# There really isnt a nicer way?
for edge in np.loadtxt('foo_1d.txt'): embedded_mark[int(edge)] = 1

#plot(f)
#plot(embedded_mark)
#interactive()

# Now let there be 1d mesh
mesh1d = EmbeddedMesh(mesh, embedded_mark, 1)

# To transfer the data
f = VertexFunction('double', mesh1d.mesh, -1)
mesh1d_indices = mesh1d.entity_map[0]

# Now this guy has data for all vertices of orig mesh
for i, value in np.loadtxt('foo_1d_data.txt'):
    # We know that this is also i-th vertex in Xd mesh and need to
    # find where it is in emebedded
    f[mesh1d_indices.index(int(i))] = value

# Let's see that we assigned correctly
for vertex in vertices(mesh1d.mesh):
    value = f[vertex]
    if value != -1:
        assert abs(vertex.point().distance(origin)-value) < 1E-13

# FIXME: EmbeddedMesh mesh should be mesh
# FIXME: at this poitn we don't have data on all the vertices as the
#        mesh is typically finer then what we started with. So we need
#        to interpolate the rest. This could be done as part of smoothing

map(os.remove, filter(lambda f: f.startswith('foo'), os.listdir('.')))
