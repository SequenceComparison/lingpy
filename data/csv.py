#! /usr/bin/env python
# *-* coding:utf-8 *-*

def loadCSV(
        infile,
        comment = '#'
        ):
    
    data = open(infile)
    out = []
    for l in data:
        if l.startswith(comment) and comment:
            pass
        else:
            out.append(
                    [k.strip() for k in l.split('\t')]
                    )
    return out

def loadDict(
        infile,
        comment='#'
        ):

    data = loadCSV(infile,comment)
    if len(data[0]) == 2:
        return dict([(l[0],l[1]) for l in data])
    else:
        return dict([(l[0],l[1:]) for l in data])

def loadList(
        infile,
        comment='#'
        ):

    data = loadCSV(infile,comment)

    return [l[0] for l in data]

class CSV(object):
    """
    
    """
    def __init__(infile):
        
        pass
