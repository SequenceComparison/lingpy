from __future__ import division
from numpy import array
from _utils import squareform

def flat_upgma(
        dict clusters,
        matrix,
        float threshold
        ):
    """
    Function recursively clusters the items in *clusters* until *threshold* is
    reached using the UPGMA algorithm (Sokal & Michener 1958).

    **Parameters:**

    * `clusters` (`dict`) -- a dictionary with integers as key and a list of
      the indices of the taxa as values
    * `matrix` -- a two-dimensional array containing the distances
    * *threshold* -- a float 
    """
    
    cdef int i,j

    # terminate when the dictionary is of length 1
    if len(clusters) == 1:
        return

    cdef list scores = []
    cdef list indices = []
    cdef list score
    cdef int valA,valB
    cdef float score_mean

    for i,valA in clusters.items():
        for j,valB in clusters.items():
            if i != j:
                score = []
                for vA in valA:
                    for vB in valB:
                        score += [matrix[vA][vB]]
                score_array = array(score)
                score_mean = score_array.mean() - 0.25 * score_array.std()
                scores.append(score_mean)
                indices.append((i,j))

    cdef float minimum = min(scores)
    
    cdef int idxA,idxB

    if minimum <= threshold:
        idxA,idxB = indices[scores.index(minimum)]
        clusters[idxA] += clusters[idxB]
        del clusters[idxB]
        return flat_upgma(
                clusters,
                matrix,
                threshold
                )
    else:
        pass

def upgma(
        dict clusters,
        matrix,
        list tree_matrix
        ):
    """
    Function recursively clusters the items in *clusters* using the UPGMA
    algorithm (Sokal & Michener 1958).
    """
    cdef int i,j,vA,vB
    cdef list valA,valB
    cdef list scores,indices,score
    cdef float score_mean

    # terminate when the dictionary is of length 1
    if len(clusters) == 1:
        return

    scores = []
    indices = []

    for i,valA in clusters.items():
        for j,valB in clusters.items():
            if i != j:
                score = []
                for vA in valA:
                    for vB in valB:
                        score += [matrix[vA][vB]]
                score_array = array(score)
                score_mean = score_array.mean() - 0.25 * score_array.std()
                scores.append(score_mean)
                indices.append((i,j))

    cdef float minimum = min(scores)
    
    cdef int idxNew = max(clusters) + 1
    
    cdef int idxA,idxB

    idxA,idxB = indices[scores.index(minimum)]
    
    clusters[idxNew] = clusters[idxA] + clusters[idxB]

    del clusters[idxA]
    del clusters[idxB]

    tree_matrix.append([idxA,idxB])
    
    return upgma(
            clusters,
            matrix,
            tree_matrix
            )

#def neighbor(
#        clusters,
#        matrix,
#        tree_matrix,
#        constant_matrix = [],
#        tracer = {}
#        ):
#    """
#    Function clusters data according to the neighbor-joining algorithm. 
#    """
#    
#    # terminate when the dictionary is of length 1
#    if len(clusters) == 1:
#        return
#    
#    # define a tracer for the order of the tree
#    if not tracer:
#        tracer = dict([(tuple([a]),b[0]) for (a,b) in clusters.items()])
#
#    # create the constant matrix when the process starts
#    if not constant_matrix:
#        constant_matrix = list(matrix.copy())
#
#    # determine the average scores
#    averages = []
#    for line in matrix:
#        averages.append(sum(line))
#    
#    # create the new matrix
#    new_matrix = matrix.copy()
#    
#    # fill in the new scores
#    for i,line in enumerate(matrix):
#        for j,score in enumerate(line):
#            if i > j:
#                new_score = (len(matrix) - 2) * score - averages[i] - averages[j]
#                new_matrix[i][j] = new_score
#                new_matrix[j][i] = new_score
#    
#    # determine the minimal score
#    scores = []
#    indices = []
#    for i in sorted(clusters.keys()):
#        for j in sorted(clusters.keys()):
#            if i < j:
#                scores.append(new_matrix[i][j])
#                indices.append((i,j))
#
#    minimum = min(scores)
#    idxA,idxB = indices[scores.index(minimum)]
#
#    # check for the average of the clusters
#    vals = []
#    for i in clusters[idxA]:
#        for j in clusters[idxB]:
#            vals.append(constant_matrix[i][j])
#    tmp_score = array(vals).mean()
#
#    idxA,idxB = indices[scores.index(minimum)]
#    
#    
#    # append the indices to the tree matrix
#    tree_matrix.append(
#            (
#                tracer[tuple(clusters[idxA])],
#                tracer[tuple(clusters[idxB])]
#                )
#            )
#
#    # create the new index for the tracer
#    idxNew = max(tracer.values()) + 1
#    tracer[tuple(clusters[idxA]+clusters[idxB])] = idxNew
#
#    # join the clusters according to the index
#    clusters[idxA] += clusters[idxB]
#    del clusters[idxB]
#
#    # create new cluster-dictionary
#    new_clusters = {}
#
#    # insert values in new dictionary
#    for i,key in enumerate(sorted(clusters.keys())):
#        new_clusters[i] = clusters[key]
#
#    # create new matrix
#    new_matrix = []
#
#    # iterate over old matrix and fill in keys for new matrix
#    for i,a in enumerate(sorted(clusters.keys())):
#        for j,b in enumerate(sorted(clusters.keys())):
#            if i < j and a != idxA and b != idxA:
#                new_matrix.append(matrix[a][b])
#            elif i < j and (a == idxA or b == idxA):
#                dist_ab = matrix[idxA][idxB]
#                dist_a = matrix[idxA][j]
#                dist_b = matrix[idxA][i]
#                new_matrix.append(((dist_a + dist_b) - dist_ab) / 2.0)
#    
#    # get values of new_clusters into clusters
#    clusters = {}
#    for key,val in new_clusters.items():
#        clusters[key] = val
#
#    # return score
#    return neighbor(
#            clusters,
#            squareform(new_matrix),
#            tree_matrix,
#            constant_matrix,
#            tracer
#            )
