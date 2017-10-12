from dolfin import Mesh
import os


# FIXME
# We don't know how to deal with tight bbox intersecting the surface and so
# are forced to have the bbox larger. Would be nice to handle this.
def mesh_around_2d(mesh_, size=1, scale=10, padding=0.05):
    '''
    From a 2d in 3d mesh (in XML format) produce a Xd mesh where
    the 2d structure is embedded. Mesh size close to strucure should 
    be size(given as multiple of hmin(), elsewhere scale * size. Padding 
    controls size of the bounding box.
    '''
    # This is a much harder problem :(
    raise NotImplementedError
