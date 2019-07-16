import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t


def calculate_distance(np.ndarray[DTYPE_t, ndim=1] v1, np.ndarray[DTYPE_t, ndim=1] v2): 

    #assert v1.dtype == DTYPE and v2.dtype == DTYPE
    
    cdef np.ndarray[DTYPE_t, ndim=1] dij
    cdef float d, a, b, c

    #dij = v1-v2
    
    a = v1[0] - v2[0]
    b = v1[1] - v2[1]
    c = v1[2] - v2[2]

    #d = np.sqrt(dij.dot(dij))

    d = (a*a + b*b + c*c)**0.5

    return d


def get_distances(np.ndarray[DTYPE_t, ndim=2] s1, np.ndarray[DTYPE_t, ndim=2] s2): 
    
    cdef list da1 = []
    cdef list da2 = []

    assert len(s1) == len(s2)
    
    #start=time.time()
    #print('Starting loop to sample distances.')

        #indexes = np.random.choice(len(sample), 2, replace=False) # <-- SO SLOW
        #indexes = [i, i + 1]
        #L = len(sample)
        #indexes = np.random.randint(0, L, size=2) # <-- MUCH FASTER (>10x)
        #while indexes[0] == indexes[1]: 
        #    indexes = np.random.randint(0, L, size=2)
    
    cdef int j, k
    cdef int lens1 = len(s1)
    cdef np.ndarray[DTYPE_t, ndim=1] v1, v2, u1, u2
    cdef float d1, d2

    for j in range(lens1): 

        #if j % 100 == 0: 
        #    print('outer loop reached event' + str(j))
            
        for k in range(j+1,lens1): 

            u1 = s1[j]
            u2 = s1[k]
            
            v1 = s2[j]
            v2 = s2[k]

            d1 = calculate_distance(u1, u2)
            d2 = calculate_distance(v1, v2)

            #distances1 = np.append(distances1, [distance1], axis=0)
            da1.append(d1)
            da2.append(d2)

    #end=time.time()
    #print('Finished getting distances. That took ' + str(end-start) + ' seconds.')
    #print('the number of distances is ' + str(len(da1)))

    return da1, da2

