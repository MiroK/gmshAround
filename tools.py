from fenics_ii.trace_tools.embedded_mesh import EmbeddedMesh
from dolfin import FunctionSpace, Function, interpolate, VertexFunction
from dolfin import dof_to_vertex_map, vertex_to_dof_map
import numpy as np


def transfer_vertex_function(mesh_fine, mesh_foo_coarse, output=VertexFunction):
    '''
    Assuming that mesh_fine is created by meshing around the mesh underlying
    mesh_foo_coarse this function interpolates the data from mesh_foo_coarse
    '''
    assert isinstance(mesh_fine, EmbeddedMesh)
    mesh_fine = mesh_fine.mesh  # FIXME: remove when EmbeddedMesh <: Mesh
    assert mesh_fine.topology().dim() == 1 and mesh_fine.geometry().dim() > 1
    
    mesh = mesh_foo_coarse.mesh()
    assert mesh.topology().dim() == 1 and mesh.geometry().dim() > 1

    # The strategy here is to interpolate into a CG1 function on mesh_fine
    # and then turn it to vertex function. NOTE: consider CG1 as function
    # for it is easier to get e.g. DG0 (midpoint values) out of it
    Vf = FunctionSpace(mesh_fine, 'CG', 1)
        
    assert mesh_foo_coarse.cpp_value_type() == 'double'
    assert mesh_foo_coarse.dim() == 0
    mesh_coarse = mesh_foo_coarse.mesh()
    Vc = FunctionSpace(mesh, 'CG', 1)
    fc = Function(Vc)
    # Fill the data
    fc.vector().set_local(mesh_foo_coarse.array()[dof_to_vertex_map(Vc)])
    fc.vector().apply('insert')

    ff = interpolate(fc, Vf)

    if output == Function: return ff
    
    # Combe back to vertex function
    vertex_foo = VertexFunction('double', mesh_fine, 0.0)
    vertex_foo.set_values(ff.vector().array()[vertex_to_dof_map(Vf)])
        
    return vertex_foo
