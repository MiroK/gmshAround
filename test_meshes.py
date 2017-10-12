from dolfin import *
from random import sample
from fenics_ii.trace_tools.embedded_mesh import EmbeddedMesh
from mesh_around_1d import mesh_around_1d
from meshconvert import convert2xml
import subprocess, os
import numpy as np


def test(path, with_mapping=True):
    '''Check the pipiline'''
    # Geo, rely on defaults
    geo, bbox = mesh_around_1d(path)
    assert os.path.exists(geo)

    gdim = 2 if bbox == 'Surface' else 3
    # Generate mesh msh
    timer = Timer('GMSH')
    timer.start()
    info('Generating %dD mesh around line. This may take some time.' % gdim)

    out =subprocess.check_output(['gmsh', '--version'], stderr=subprocess.STDOUT)
    assert out.split('.')[0] == '3', 'Gmsh 3+ is required'
    
    ccall= 'gmsh -%d -optimize %s' % (gdim, geo)
    subprocess.call(ccall, shell=True)
    
    root, ext = os.path.splitext(geo)
    msh_file = '.'.join([root, 'msh'])
    assert os.path.exists(msh_file)
    info('Done in %g s' % timer.stop())
    
    # Convert to xml
    xml_file = '.'.join([root, 'xml'])
    convert2xml(msh_file, xml_file)

    # There should be a mesh, physical, embedded
    assert os.path.exists(xml_file)
    volume_file = '.'.join(['_'.join([root, 'physical_region']), 'xml'])
    assert os.path.exists(volume_file)

    delete = [xml_file, volume_file]
    if gdim == 2:
        manifold_file = '.'.join(['_'.join([root, 'facet_region']), 'xml'])
    else:
        manifold_file = '.'.join(['_'.join([root, 'edge_region']), 'xml'])
    assert os.path.exists(manifold_file)
    delete.append(manifold_file)

    # Check the consistency, vertex - vertex correspondence. Gmsh kept
    # the points in 3d as we were feeding it on 1d. So we can find vertices
    # in 3d mesh no problem and EmbeddedMesh has a map from its vertices
    # to 3d vertices.
    meshXd = Mesh(xml_file)
    mesh1d = Mesh(path)

    assert np.linalg.norm(mesh1d.coordinates() -
                          meshXd.coordinates()[:mesh1d.num_vertices()]) < 1E-13
    
    # Cleanup
    delete.extend([geo, msh_file])
    map(os.remove, delete)

    return True


def mesh_2d():
    '''Generate data'''
    mesh = UnitSquareMesh(32, 32)
    mesh.init(1)

    n = mesh.num_entities(1)
    f = FacetFunction('size_t', mesh, 0)

    for i in sample(range(n), n/4): f[i] = 1
    
    mesh = EmbeddedMesh(mesh, f, 1).mesh

    out = 'mesh2d.xml'
    File(out) << mesh

    return out


def mesh_3d():
    '''Generate data'''
    mesh = UnitCubeMesh(16, 16, 16)
    mesh.init(1)

    n = mesh.num_entities(1)
    f = EdgeFunction('size_t', mesh, 0)

    for i in sample(range(n), n/4): f[i] = 1

    mesh = EmbeddedMesh(mesh, f, 1).mesh

    out = 'mesh3d.xml'
    File(out) << mesh

    return out

# --------------------------------------------------------------------

if __name__ == '__main__':
    # Markers -> geo -> gmsh -> assumptions -> cleanup
    assert test(mesh_2d())
    assert test(mesh_3d())
