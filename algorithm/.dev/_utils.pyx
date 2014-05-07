from __future__ import division
from numpy import zeros

cdef extern from "math.h":
    float sqrt(float theta)

cdef class LingpyArray(object):
    """
    An extension of the numpy array object which allows the storage of lists in
    two-dimensional arrays.
    """

    def __cinit__(
            self,
            list input_list
            ):
        
        cdef int h = len(input_list)
        cdef int w = len(input_list[0])

        self.array = zeros((h,w),dtype='int')

        self.dictionary = {}

        cdef int count = 0
        for i,line in enumerate(input_list):
            for j,itm in enumerate(line):
                self.dictionary[count] = input_list[i][j]
                self.array[i][j] = count
                count += 1

    def __getitem__(self,x):

        ind = self.array[x]

        try:
            return self.dictionary[self.array[x]]
        except:
            try:
                return [self.dictionary[i] for i in ind]
            except:
                return [[self.dictionary[i] for i in j] for j in ind]

def squareform(list x):
    """
    A simplified version of the squareform function in scipy.spatial.distance.
    The version is less safe than the original one, yet it is much
    faster.
    """
    cdef int i,j,k
    cdef int l = len(x)

    # calculate the length of the square
    cdef int s = int(sqrt(2 * l) + 1)
    
    out = zeros((s,s))
    
    k = 0
    for i in range(s):
        for j in range(s):
            if i < j:
                out[i][j] = x[k]
                out[j][i] = x[k]
                k += 1
    return out
