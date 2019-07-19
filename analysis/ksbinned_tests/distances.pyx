# cython: profile=True
# cython: boundscheck=False, wraparound=False, nonecheck=False


import numpy as np
cimport numpy as np
from libc.math cimport exp
from libc.math cimport log


DTYPE = np.double
ctypedef np.double_t DTYPE_t


def calculate_distance(double v1, double v2, double v3, double u1, double u2, double u3): 

    #assert v1.dtype == DTYPE and v2.dtype == DTYPE
    
    #cdef np.ndarray[DTYPE_t, ndim=1] dij
    cdef double d, a, b, c, psi

    #dij = v1-v2
    
    a = v1 - u1
    b = v2 - u2
    c = v3 - u3

    #d = np.sqrt(dij.dot(dij))

    d = (a*a + b*b + c*c)

    psi = exp(- d / 0.6**2)

    return psi


def get_distances(np.ndarray[DTYPE_t, ndim=2] s1, np.ndarray[DTYPE_t, ndim=2] s2): 
    
    cdef int lens1 = len(s1)

    cdef list da1 = [None]*((lens1*lens1 - lens1)/2)
    cdef list da2 = [None]*((lens1*lens1 - lens1)/2)

    cdef int j, k, i
        #cdef np.ndarray[DTYPE_t, ndim=1] v1, v2, u1, u2
    cdef double d1, d2
    cdef double u1, u2, u3, v1, v2, v3, w1, w2, w3, x1, x2, x3
    
    i = 0
    for j in range(lens1): 

        #if j % 100 == 0: 
        #    print('outer loop reached event' + str(j))
        u1 = s1[j][0]
        u2 = s1[j][1]
        u3 = s1[j][2]

        w1 = s2[j][0]
        w2 = s2[j][1]
        w3 = s2[j][2] 

            
        for k in range(j+1,lens1): 
            
            v1 = s1[k][0]
            v2 = s1[k][1]
            v3 = s1[k][2]

            x1 = s2[k][0]
            x2 = s2[k][1]
            x3 = s2[k][2]

            d1 = calculate_distance(u1, u2, u3, v1, v2, v3)
            d2 = calculate_distance(w1, w2, w3, x1, x2, x3)

            #distances1 = np.append(distances1, [distance1], axis=0)
            #da1.append(d1)
            #da2.append(d2)
            da1[i] = d1
            da2[i] = d2

            i=i+1

    #end=time.time()
    #print('Finished getting distances. That took ' + str(end-start) + ' seconds.')
    #print('the number of distances is ' + str(len(da1)))

    return da1, da2

