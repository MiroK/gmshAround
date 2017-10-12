from dolfin import Mesh, cells, edges
import os

# FIXME
# We don't know how to deal with tight bbox intersecting the surface and so
# are forced to have the bbox larger. Would be nice to handle this.

def mesh_around_1d(mesh, size=1, scale=10, padding=0.05):
    '''
    From a 1d in xd (X > 1) mesh (in XML format) produce a Xd mesh where
    the 1d structure is embedded. Mesh size close to strucure should 
    be size(given as multiple of hmin(), elsewhere scale * size. Padding 
    controls size of the bounding box.
    '''
    dot = mesh.find('.')
    root, ext = mesh[:dot], mesh[dot:]
    assert ext == '.xml' or ext == '.xml.gz', ext

    mesh = Mesh(mesh)
    gdim = mesh.geometry().dim()
    assert gdim > 1 and mesh.topology().dim() == 1

    x = mesh.coordinates()
    mesh.init(1, 0)

    # Compute fall back mesh size:
    assert size > 0
    size = mesh.hmin()*size

    # Don't allow zero padding - collision of lines with bdry segfaults
    # too ofter so we prevent it
    assert padding > 0
    # Finally scale better be positive
    assert scale > 0

    point = (lambda xi: tuple(xi) + (0, ))\
            if gdim == 2 else (lambda xi: tuple(xi))
    
    geo = '.'.join([root, 'geo'])
    with open(geo, 'w') as outfile:
        # Setup
        outfile.write('SetFactory("OpenCASCADE");\n')
        outfile.write('size = %g;\n' % size)
        outfile.write('SIZE = %g;\n' % (size*scale))

        # Points
        fmt = 'Point(%d) = {%.16f, %.16f, %.16f, size};\n'
        for i, xi in enumerate(x, 1):
            outfile.write(fmt % ((i, ) + point(xi)))
        # Lines
        fmt = 'Line(%d) = {%d, %d};\n'
        for i, cell in enumerate(cells(mesh), 1):
            outfile.write(fmt % ((i, ) + tuple(cell.entities(0)+1)))

        # BBox
        xmin, xmax = x.min(0), x.max(0)
        padding = (xmax-xmin)*padding/2.
        xmin -= padding
        xmax += padding
        dx = xmax - xmin

        if gdim == 2 or dx[-1] < 1E-14: # All points are on a plane
            rect = 'Rectangle(1) = {%g, %g, %g, %g, %g};\n' % (xmin[0],
                                                               xmin[1],
                                                               0 if gdim == 2 else xmin[2],
                                                               dx[0],
                                                               dx[1])
            outfile.write(rect)
            bbox = 'Surface'
        else:
            box = 'Box(1) = {%g, %g, %g, %g, %g, %g};\n' % (xmin[0],
                                                            xmin[1],
                                                            xmin[2],
                                                            dx[0],
                                                            dx[1],
                                                            dx[2])
            outfile.write(box)
            bbox = 'Volume'
        
        # Crack
        for line in xrange(1, mesh.num_cells()+1):
            outfile.write('Line{%d} In %s{1};\n' % (line, bbox))
            
        # Add Physical volume/surface
        outfile.write('Physical %s(1) = {1};\n' % bbox)

        # Add Physical surface/line
        lines = ', '.join(map(lambda v: '%d' % v, xrange(1, mesh.num_cells()+1)))
        outfile.write('Physical Line(1) = {%s};\n' % lines)
    return geo, gdim
