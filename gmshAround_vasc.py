# This script illustrates how the Xd-1d mesh could be makde from 1d
from dolfin import Mesh, MeshFunction, Timer
from mesh_around_1d import mesh_around_1d
from meshconvert import convert2xml
import subprocess, os
import numpy as np

timer = Timer('Meshing')
timer.start()

path = 'vasc_mesh.xml.gz'
geo, _ = mesh_around_1d(path, size=1.)

data = MeshFunction('int', Mesh(path), 'widths.xml.gz')

out =subprocess.check_output(['gmsh', '--version'], stderr=subprocess.STDOUT)
assert out.split('.')[0] == '3', 'Gmsh 3+ is required'

ccall= 'gmsh -3 -optimize %s' % geo
subprocess.call(ccall, shell=True)
    
# Convert
xml_file = 'vasc_GMSH.xml'
msh_file = 'vacs_GMSH.msh'
convert2xml(msh_file, xml_file)

# Throw away the line function as it is expensive to store all the edges
meshXd = Mesh(xml_file)
# Make sure that the assumption of correspondence hold
assert np.linalg.norm(mesh.coordinates() -
                      meshXd.coordinates()[:mesh.num_vertices()]) < 1E-13

line_f = MeshFunction('size_t', meshXd, 'vasc_GMSH_edge_region.xml')
    
# Just a list of edges which are 1
tag_file = 'vasc_GMSH_vessel_tags.txt'
np.savetxt(tag_file, [edge.index() for edge in SubsetIterator(line_f, 1)])

# Transfer is not done here because I don't wantt to create embedded
# mesh etc
dt = timer.stop()
print 'Done in %g s' % dt, 'look for', xml_file, 'and', tag_file
