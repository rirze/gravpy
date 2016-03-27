#!/usr/bin/python

import numpy as np
import numpy.linalg as la

def sub(u,v):
    '''Subtract two vectors '''
    return (u[0]-v[0],u[1]-v[1])

def det(u,v):
    '''Return the determinant of the 2-D matrix formed by u and v'''
    return u[0]*v[1]-u[1]*v[0]

def inside_triangle(P,tri):
    '''Explicit expression of the components in the algebra of finding barycentric coordiantes. Plugging in an x and y value for v0,v1,v2 reduces the computation to three determinants from before which computed six dot products.'''
    # Code derived from http://www.blackpawn.com/texts/pointinpoly/
    # Uses barycentric coordinates to find whether the point is inside the triangle

    A,B,C = tri
    
    v0 = sub(C,A)
    v1 = sub(B,A)
    v2 = sub(P,A)

    invDenom = 1.0/det(v1,v0)
    
    u = det(v2,v1) * (-invDenom)
    v = det(v2,v0) *   invDenom

    return u >= 0 and v >= 0 and u + v < 1

def find(v,triangles):
    '''Given a point v and a list of triangles defined by three verticies, this function will return the indicies of the triangles in which the point is inside the triangle.'''

    indices = [i for i,triangle in enumerate(triangles) if inside_triangle(v,triangle)]
            
    return indices

def custom_2D_det(u,v):
    u1,u2 = u.T
    v1,v2 = v.T
    return u1*v2-u2*v1

    
def np_inside_triangle(A,B,C,P):
    '''Explicit expression of the components in the algebra of finding barycentric coordiantes. Plugging in an x and y value for v0,v1,v2 reduces the computation to three determinants from before which computed six dot products. This particular function uses numpy vector operations, where the inputs are vectors of pairs. 'a,b,c' are the lists of pairs of the triangle vertices, and p is still a single point for which we are searching for.'''
    # Code derived from http://www.blackpawn.com/texts/pointinpoly/
    # Uses barycentric coordinates to find whether the point is inside the triangle
    # Uses numpy vector operations

    v0 = C-A
    v1 = B-A
    v2 = P-A
    
    with np.errstate(divide='ignore',invalid='ignore'): 
        #my understanding is that invalid values arise from floating point errors, so nothing to worry about? still gives good results...
        invDenom = 1.0/custom_2D_det(v1,v0) #divide/0 
    
        u = custom_2D_det(v2,v1) * (-invDenom) #invalid value in det
        v = custom_2D_det(v2,v0) *   invDenom  #invalid value in det
        
        coords = np.vstack((u,v,1-u-v))
        indices = np.nonzero(np.all(coords>=0,axis=0)) #invalid value in greater_equal
    
    return indices

def find2(p,triangles):
    '''Wrapper for the algorithm that uses numpy vector operations. Use this one for better speed!'''
    a,b,c = np.transpose(triangles,[1,0,2])
    result = np_inside_triangle(a,b,c,p)
    
    if result[0].size == 0:
        raise ValueError("Point {} not found in given triangles".format(p))
    else:
        return result

#****Speedups****
# dot products -> determinants in inside_triangle(): 80 msec to 60 msec
# for loop -> list comprehension in find(): 60 msec to 56 msec
# find() -> find2() -- full transition to numpy vector operations: 56 msec to 15 msec
# np.det() -> custom_2D_det(): 15 msec to 1 msec 

    
