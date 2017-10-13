from dolfin import Mesh, MeshEditor, cells
import networkx as nx


def assertions_hold(mesh):
    '''Line mesh in R^d, d > 1'''
    return mesh.topology().dim() == 1 and mesh.geometry().dim() > 1


def branch_terminals(mesh):
    '''
    In our definition a branch connects two vertices of the mesh which
    are such that either a single cell is connected to the vertex or the
    vertex is a bifurcation = 3 and more cells share it. Here we return
    a list (index) of such points and the cells that them
    '''
    assert assertions_hold(mesh)

    mesh.init(0, 1)
    mesh.init(1, 0)
    v2c = mesh.topology()(0, 1)

    mapping = dict()
    for v in xrange(mesh.num_vertices()):
        v_cells = v2c(v)
        if len(v_cells) == 1 or len(v_cells) > 2:
            mapping[v] = set(v_cells.tolist())

    return mapping


def find_branches(mesh):
    '''
    Produces a cell function marking each branch of the mesh with different
    color.
    '''
    terminals_map = branch_terminals(mesh)

    v2c = mesh.topology()(0, 1)
    c2v = mesh.topology()(1, 0)

    branches, terminals = [], []
    while terminals_map:
        start = vertex = next(iter(terminals_map))
        # Visited vertices of the branch
        visited = set((start, ))
        edge = terminals_map[start].pop()
        # Edges in the branch
        branch = [edge]
        while True:
            v0, v1 = c2v(edge)
            vertex = v0 if v1 in visited else v1
            visited.add(vertex)
            # Terminal vertex; remove how I got here
            if vertex in terminals_map:
                # How I got here
                terminals_map[vertex].remove(edge)
                if not terminals_map[vertex]: del terminals_map[vertex]
                if not terminals_map[start]: del terminals_map[start]
                break
            # Otherwise, compute new edge
            try:
                e0, e1 = v2c(vertex)
                edge = e0 if e1 in branch else e1
            except ValueError:
                edge = v2c(vertex)[0]

            assert edge not in branch
            branch.append(edge)
            
        branches.append(branch)
        terminals.append((start, vertex))

    return terminals, branches


def reduced_mesh(mesh):
    '''
    Represent each branch only as a segment / single cell in the mesh.
    The mesh has a map from new mesh to old mesh vertices
    '''
    terminals, _ = find_branches(mesh)

    # Unique nodes, this is map from rmesh nodes to parent
    nodes = map(int, set(sum(terminals, ())))
    # Let's make reverse for purpose of creating the mesh
    old2new = {old: new for new, old in enumerate(nodes)}

    # Their coordinates
    x = mesh.coordinates()[nodes]

    rmesh = Mesh()
    editor = MeshEditor()

    editor.open(rmesh, 1, mesh.geometry().dim())
    editor.init_vertices(len(x))
    editor.init_cells(len(terminals))

    # Add vertices
    for vi, v in enumerate(x): editor.add_vertex(vi, v)

    # Add cells
    for ci, c in enumerate(terminals):
        editor.add_cell(ci, *map(lambda v: old2new[v], c))

    editor.close()

    # How to do this with MeshData
    rmesh.parent_vertex_indices = nodes
    
    return rmesh


class Walker(object):
    def __init__(self, mesh):
        assert assertions_hold(mesh)

        self.v2c = mesh.topology()(0, 1)
        self.c2v = mesh.topology()(1, 0)
        self.volumes = [cell.volume() for cell in cells(mesh)]

    def walk(self, branch, endpoints=None):

        v2c = self.v2c
        c2v = self.c2v

        if endpoints is None:
            # The one that is not connected to remaining cells
            cell = branch[0]
            f0, f1 = c2v(cell)
            first = f1 if any(c in branch for c in set(v2c(f0))-set((cell, ))) else f0

            cell = branch[-1]
            f0, f1 = c2v(cell)
            last = f1 if any(c in branch for c in set(v2c(f0))-set((cell, ))) else f0

            endpoints = [first, last]

        first, last = endpoints
        # Consistency
        assert first in c2v(branch[0]) and last in c2v(branch[-1])
        
        # Is this really a branch
        assert len(v2c(first)) == 1 or len(v2c(first)) > 2 
        assert len(v2c(last)) == 1 or len(v2c(last)) > 2

        yield first, 0

        volumes = self.volumes
        
        d = 0
        prev = first
        for cell in branch:
            v0, v1 = c2v(cell)
            link = v0 if v1 == prev else v1

            d += volumes[cell]
            
            yield link, d

            prev = link
            
# -------------------------------------------------------------------

from dolfin import MeshFunction, File, CellFunction, Cell, near

mesh = Mesh('vasc_mesh.xml.gz')
data = MeshFunction('int', mesh, 'widths.xml.gz')
data = data.array()

# Original
mesh.init(1, 0)
data_o = CellFunction('double', mesh, 0)

for cell in cells(mesh):
    data_o[cell] = sum(data[cell.entities(0)])/2.
File('data_original.pvd') << data_o

terminals, branches = find_branches(mesh)

rmesh = reduced_mesh(mesh)
node_map = rmesh.parent_vertex_indices

rmesh.init(1, 0)
# Let's see how bumpy the coarse representation is
rmesh_foo = CellFunction('double', rmesh, 0)
for cell in cells(rmesh): rmesh_foo[cell] = sum(data[cell.entities(0)])/2.

File('data_coarse.pvd') << rmesh_foo

w = Walker(mesh)

# Checks: distance increases monot, all the points produced by the walker
# are all the points on the branch
if False:
    mesh.init(1, 0)
    for branch in branches:

        pts = set()
        for p, d in w.walk(branch):
            pts.add(p)
        distance = d

        distance_ = sum(Cell(mesh, c).volume() for c in branch)
        assert near(distance, distance_), (distance, distance_)

        pts_ = set(sum((Cell(mesh, c).entities(0).tolist() for c in branch), []))
        assert pts ==  pts_, (pts, pts_)

import matplotlib.pyplot as plt
# Let's see how badly the data oscillates
for branch in branches:
    x, y = [], []
    for p, d in w.walk(branch):
        x.append(d)
        y.append(data[p])
    plt.plot(x, y)
plt.show()
