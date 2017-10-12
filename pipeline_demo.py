# This script illustrates how the Xd-1d mesh could be makde from 1d
from dolfin import *
from random import sample
from fenics_ii.trace_tools.embedded_mesh import EmbeddedMesh
from mesh_around_1d import mesh_around_1d
from tools import transfer_vertex_function
from meshconvert import convert2xml
import subprocess, os
import numpy as np

dim = 3
if dim == 2:
    # Fabricate some 1d mesh and store it together with some vertex function
    mesh = UnitSquareMesh(32, 32)
    mesh.init(1)

    n = mesh.num_entities(1)
    f = FacetFunction('size_t', mesh, 0)

    for i in sample(range(n), n/4): f[i] = 1
else:
    assert dim == 3

    mesh = UnitCubeMesh(16, 16, 16)
    mesh.init(1)

    n = mesh.num_entities(1)
    f = EdgeFunction('size_t', mesh, 0)

    for i in sample(range(n), n/8): f[i] = 1
        
# The 1d mesh
mesh = EmbeddedMesh(mesh, f, 1).mesh
# Fake data
data = VertexFunction('double', mesh, 0)

# Set to something linear depending on coords i) can check ii) can interpolate
# exactly
origin = Point(0., 0.)
for vertex in vertices(mesh):
    data[vertex] = sum(vertex.point().array())

out, out_data = 'foo.xml', 'foo_data.xml'
File(out) << mesh
File(out_data) << data

# Create the 2d mesh around it
geo, _ = mesh_around_1d(out)

out =subprocess.check_output(['gmsh', '--version'], stderr=subprocess.STDOUT)
assert out.split('.')[0] == '3', 'Gmsh 3+ is required'

ccall= 'gmsh -%d -optimize %s' % (dim, geo)
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

if dim == 2:
    line_f = MeshFunction('size_t', meshXd, 'foo_facet_region.xml')
else:
    line_f = MeshFunction('size_t', meshXd, 'foo_edge_region.xml')
    
assert all(line_f[e] == 1 or line_f[e] == 0 for e in edges(meshXd))

# Just a list of edges which are 1
np.savetxt('foo_1d.txt', [edge.index() for edge in SubsetIterator(line_f, 1)])

################
# Loading part
################
# Now we going to build the 1d mesh back and get the data as well
mesh = Mesh('foo.xml')  # 2d ro 3d

# FIXME: embedded mesh should accept a list
embedded_mark = EdgeFunction('size_t', mesh, 0)
# There really isnt a nicer way?
for edge in np.loadtxt('foo_1d.txt'): embedded_mark[int(edge)] = 1

#plot(f)
#plot(embedded_mark)
#interactive()

# Now let there be 1d mesh
mesh1d = EmbeddedMesh(mesh, embedded_mark, 1)

# Transfer the data from XML on original 1d to here
f = transfer_vertex_function(mesh1d, data)

# Let's we that we did it okay
for vertex in vertices(mesh1d.mesh):
    assert abs(sum(vertex.point().array()) - f[vertex]) < 1E-13, \
        ('%.16f' % abs(f[vertex] - sum(vertex.point().array())))
    
map(os.remove, filter(lambda f: f.startswith('foo'), os.listdir('.')))

# FIXME: EmbeddedMesh mesh should be mesh
# FIXME: smoothing
