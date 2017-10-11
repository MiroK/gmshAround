from dolfin import Mesh, cells
import os


def mesh_around_1d(mesh, size=-1, scale=10, padding=0.05):
    '''
    From a Xd in 1d (X > 1) mesh (in XML format) produce a Xd mesh where
    the 1d structure is embedded. Mesh size close to strucure should 
    be size, elsewhere scale * size. Padding controls size of the bounding
    box.
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
    if size < 0:
        size = mesh.hmin()/2.
    else:
        assert size > 0

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
    return geo, bbox

# -------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()
    parser.add_argument('mesh', type=str, help='path to XML mesh')
    parser.add_argument('--size', type=float, default=-1,
                        help='characteristic mesh size at 1d points (default is h/2)')
    parser.add_argument('--scale', type=float, default=10.,
                        help='scale * size will be the charactic size elsewhere')
    parser.add_argument('--padding', type=float, default=0.05,
                        help='multiple of bbox size to to bbox (default 0.1)')
    
    args = parser.parse_args()

    sys.exit(mesh_around_1d(args.mesh, args.size, args.scale, args.padding))
