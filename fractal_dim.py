from math import log
import numpy as np
import dolfin as df
import operator


def estimate_from_convergence(y, x):
    '''Half step estimate'''
    assert len(x) == len(y)
    
    if len(y) >= 2:
        return -log(y[-1]/float(y[-2]))/log(x[-1]/x[-2])
    else:
        return np.nan

    
def least_square_estimate(y, x):
    '''Fit for y = x^(-p)+a'''
    assert len(x) == len(y)

    if len(y) >= 2:
        A = np.vstack([np.log(x), np.ones(len(x))]).T
        p, _ = np.linalg.lstsq(A, np.log(y))[0]
        return -p
    else:
        return np.nan

    
def collides(box, mesh):
    b = box
    for cell in tree:
        # Collision with b
    pass


def split(box):
    '''Split the box into 8 children'''
    x, y = box
    m = (x+y)/2
    return [(np.array([x[0], x[1], x[2]]), np.array([m[0], m[1], m[2]])),
            (np.array([x[0], m[1], x[2]]), np.array([m[0], y[1], m[2]])),
            (np.array([m[0], x[1], x[2]]), np.array([y[0], m[1], m[2]])),
            (np.array([m[0], m[1], x[2]]), np.array([y[0], y[1], m[2]])),
            (np.array([x[0], x[1], m[2]]), np.array([m[0], m[1], y[2]])),
            (np.array([x[0], m[1], m[2]]), np.array([m[0], y[1], y[2]])),
            (np.array([m[0], x[1], m[2]]), np.array([y[0], m[1], y[2]])),
            (np.array([m[0], m[1], m[2]]), np.array([y[0], y[1], y[2]]))]


def fractal_dim(mesh):
    '''
    Produce estimates of fractal dimension of the mesh by box counting
    '''
    assert mesh.topology().dim() == 1 and mesh.geometry().dim() == 3

    root = bbox(mesh)
    
    gdim = mesh.geometry().dim()
    sizes = root[1] - root[0]
    size = reduce(operator.__mult__, sizes)**1./3

    leafs = [root]
    while True:
        count = 0
        for box in leafs:
            if collides(box, cells(mesh)):
                count += 1
                new_leafs.extend(split(box))
                
        yield count, size

        leafs = new_leafs
        size /= 2.
        
# -------------------------------------------------------------------

#if __name__ == '__main__':


    
    
mesh = df.Mesh('vasc_mesh.xml.gz')

N_history, eps_history = [], []

print bbox(mesh)

print bbox(next(df.cells(mesh)))

# i = 0
# for N, eps in fractal_dim(mesh):
#     N_history.append(N)
#     eps_history.append(eps)

#     print i
#     print N, eps
#     print '\t current', -log(N)/log(eps)
#     print '\t cvrg', estimate_from_convergence(N_history, eps_history)
#     print '\t lstsq', least_square_estimate(N_history, eps_history)
#     print        
