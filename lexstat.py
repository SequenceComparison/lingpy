# *-* coding: utf-8 *-*
"""
This is the basic module for the 
"""

from __future__ import division,print_function
import os
from numpy import array,zeros,log2 
from pickle import load,dump
from data import *
from algorithm import *
from algorithm.cluster import _flat_upgma,_flat_neighbor,_upgma,_neighbor
from align.multiple import _Multiple
from output.color import colorRange

class LexStat(object):
    """
    Basic class for handling lexicostatistical datasets.

    Parameters
    ----------
    infile : file
        A file in ``lxs``-format.  
    
    Notes
    -----
    The LexStat class serves as the base class for the handling of
    lexicostatistical datasets (see :evobib:`Swadesh1955` for a detailed description of the
    method of lexicostatistics). It provides methods for data conversion, when
    analyses on cognacy have been conducted in a qualitative way, and also
    allows to carry out cognate judgments automatically, based on the different
    methods described in :evobib:`List2012a`.

    The input data for LexStat is a simple tab-delimited text file with the
    language names in the first row, an ID in the first column, and the data in
    the columns corresponding to the language names. Additionally, the file can
    contain headwords corresponding to the IDs and cognate-IDs, specifying
    which words in the data are thought to be cognate.  This structure is
    almost the same as the one employed in the Starling database program (see
    http://starling.rinet.ru). Synonyms are also specified in the same way by
    simply adding additional rows with the same ID. The following is an example
    for the possible structure of an input file::
    
        ID  Word    German  COG   English   COG ...
        1   hand    hantʰ   1     hæːnd     1   ...
        2   fist    faustʰ  2     fist      2   ...
        ... ...     ...     ...   ...       ... ...  

    .. todo:: csv-documentation
       
       Add a description of the current csv-documentation for interaction with
       qlc-library.
    """
    def __init__(
            self,
            infile,
            ignore_errors = True,
            verbose = False
            ):
        
        # switch of printing etc.
        self.verbose = verbose

        # store the name of the input file
        if infile.endswith('.lxs'):
            self.infile = os.path.split(infile)[1].replace('.lxs','')
            self._init_lxs(infile,ignore_errors)
        elif infile.endswith('.csv'):
            self.infile = os.path.split(infile)[1].replace('.csv','')
            self._init_csv(infile,ignore_errors)
        else:
            self.infile = os.path.split(infile)[1]
            self._init_csv(infile,ignore_errors)

    def _init_csv(
            self,
            infile,
            ignore_errors = True,
            ):
        """
        Basic routine for input files in csv-format.
        """
        # load the data
        try:
            txt = array(loadtxt(infile),dtype='str')
        except IOError:
            txt = array(loadtxt(infile+'.csv'),dtype='str')

        # carry out primary check of data, such as, e.g. the search for
        # multiple spaces
        for i,line in enumerate(txt):
            if len([0 for j in line if '  ' in j]) > 0:
                if self.verbose: print("[!] Warning, data contains multiple spaces in line {0}!".format(i+1)
                        )

        # converter for headlines into internal format
        converter = {
                "ipa" : "words",
                "gloss" : "items",
                "item" : "items",
                "prostrings" : "prostrings",
                "taxon" : "taxa",
                "language" : "taxa",
                "languages" : "taxa",
                "iso" : "taxa",
                "taxa" : "taxa",
                "numbers" : "numbers",
                "tokens" : "tokens",
                "glossid" : "numbers",
                "gloss_id" : "numbers",
                "id" : "ids",
                "classes" : "classes",
                "cogs" : "cogs",
                "cognates" : "cogs",
                "cognate_id" : "cogs",
                "cog_id" : "cogs",
                "cogid" : "cogs"
                }
        # store the index of the columns in the structure variable of the input
        # file
        self._ifs = {}
        for i,x in enumerate(txt[0]):
            y = x.strip()
            try:
                self._ifs[converter[y.lower()]] = (i,x)
            except KeyError:
                if y.lower() not in converter.values():
                    self._ifs[y.lower()] = (i,x)
                else:
                    self._ifs[y.lower()+'-x'] = (i,x)

        # check infile for required information
        required = [
                ("words","ipa"),
                ("numbers","gloss_id"),
                ("taxa","language"),
                ("ids","ID")
                ]

        for r in required:
            if r[0] not in self._ifs.keys():
                if self.verbose: print('[!] Missing "{0}"-column in input file!'.format(r[1]))
        
        # store all values of txt in an array
        data = array(txt[1:])
        
        # get the taxa XXX note that taxa are now sorted XXX
        self.taxa = sorted(list(
                set(
                    [
                        x.strip().strip('.') for x in data[:,self._ifs['taxa'][0]]
                        ]
                        )
                    ))
        self.width = len(self.taxa)
        
        ## get the gloss_ids
        #self.numbers = [int(x) for x in data[:,self._ifs['numbers'][0]]]

        ## create the items, if they are not in the data
        #if 'items' not in self._ifs:
        #    self.items = [str(x) for x in self.numbers]
        #else:
        #    self.items = data[:,self._ifs['items'][0]]

        # create idxd first, i.e. a 2d-array in which each id corresponds to the
        # later position in the two-dimensional array
        
        # create a temporary dictionary that groups all language entries according to identical
        # meanings, and a temporary dictionary which groups all lines by ID
        tmp,dataById = {},{}
        
        for line in data:
            lng = self.taxa.index(
                    line[self._ifs['taxa'][0]].strip().strip('.')
                    )
            idx = int(line[self._ifs['ids'][0]])
            wrd = line[self._ifs['words'][0]].strip()
            gls = int(line[self._ifs['numbers'][0]])

            try:
                tmp[gls] += [(idx,lng,wrd)]
            except KeyError:
                tmp[gls] = [(idx,lng,wrd)]

            dataById[idx] = line
        
        # create idxd, the word-array and the array for the ids
        self.words = [] # array for the words
        self.idxs = [] # array for the ids
        self.idxd = {}
        self.idx = {}
        self.dataById = dataById # input data in straight format
        
        # add first line to dataById
        self.dataById[0] = txt[0]
        
        lines = 0
        for key,vals in sorted(tmp.items(),key=lambda x: x[0]):
            # get the maximum number of synonyms in the data
            m = max([[i[1] for i in vals].count(j) for j in range(self.width)])
            tmp_words= [['-' for i in range(self.width)] for j in range(m)]
            tmp_idxs = [[0 for i in range(self.width)] for j in range(m)]
            
            # fill in the values for each line
            for val in vals:
                k = 0
                finished = False
                while not finished:
                    if tmp_idxs[k][val[1]]:
                        k += 1
                        
                    else:
                        tmp_idxs[k][val[1]] = val[0]
                        tmp_words[k][val[1]] = val[2]
                        self.idxd[val[0]] = (lines+k,val[1])
                        finished = True
    
            # fill in values for idx
            self.idx[key] = [i+lines for i in range(m)]

            # increase the line-counter
            lines += m
            
            # fill in the values for the words and the ids
            for w in tmp_words:
                self.words.append(w)
            for i in tmp_idxs:
                self.idxs.append(i)

        # convert words and idxs to arrays
        self.words = array(self.words)
        self.idxs = array(self.idxs)
        
        # create self.numbers
        self.numbers = []
        self.items = []

        for i,line in enumerate(self.words):
            tmp_numbers = []
            tmp_items = []
            for j,cell in enumerate(line):
                if cell != '-':
                    tmp_numbers.append(
                            self.dataById[self.idxs[i][j]][
                                self._ifs['numbers'][0]
                                ]
                            )
                    if 'items' in self._ifs:
                        tmp_items.append(
                                self.dataById[self.idxs[i][j]][
                                    self._ifs['items'][0]
                                    ]
                                )
                    else:
                        tmp_items.append(str(self.numbers[-1]))
            self.numbers.append(tmp_numbers[0])
            self.items.append(tmp_items[0])


        # create self.cogs, if this is not already there
        if 'cogs' not in self._ifs:
            self.cogs = zeros(self.words.shape,dtype='int')
        else:
            self.cogs = zeros(self.words.shape,dtype='int')
            for i,line in enumerate(self.cogs):
                for j,cell in enumerate(line):
                    if self.idxs[i][j] != 0:
                        self.cogs[i][j] = int(
                                self.dataById[self.idxs[i][j]][self._ifs['cogs'][0]]
                                )

        # create a specific format string in order to receive taxa of equal
        # length
        mtax = max([len(t) for t in self.taxa])
        self._txf = '{0:.<'+str(mtax)+'}'

        # get the prosodic strings and the context information
        errors = {'chars':{},'words':{}} # list to store all errors in the document
        self.tokens = self.words.copy().tolist()
        self.prostrings = self.words.copy()
        self.alignments = self.words.copy().tolist()

        # we add the sonority profiles since they are needed for
        # alignment-based clustering XXX
        self.sonars = self.words.copy().tolist()

        for i in range(self.words.shape[0]):
            for j in range(self.words.shape[1]):
                if self.words[i][j] != '-':
                    try:
                        if 'tokens' in self._ifs:
                            tokens = unicode(
                                    dataById[self.idxs[i][j]][self._ifs['tokens'][0]],
                                    'utf-8'
                                    ).split(' ')
                        else:
                            tokens = ipa2tokens(self.words[i][j])
                        sonar = [int(x) for x in tokens2class(tokens,art)]
                        prostring = prosodic_string(sonar)
                        self.tokens[i][j] = tokens
                        self.prostrings[i][j] = prostring

                        # XXX be careful with list assignment and remember always
                        # XXX to copy lists properly in python!
                        self.alignments[i][j] = [t for t in tokens] 
                        self.sonars[i][j] = sonar
                    except:
                        errstring = "[!] Error with '{0}' in <{1}> (ID: {2}!"
                        if tokens:
                            check = check_tokens(tokens)
                            if check:
                                checkstring = ' '.join([err[1] for err in check])
                                if not ignore_errors:
                                    if self.verbose: print(errstring.format(
                                        checkstring.encode('utf-8'),
                                        self.words[i][j],
                                        self.idxs[i][j]
                                        ))
                                for idx,err in check:
                                    try:
                                        errors['chars'][err] += 1
                                    except KeyError:
                                        errors['chars'][err] = 1
                            else:
                                checkstring = ' '.join(tokens).encode('utf-8')
                                if not ignore_errors:
                                    if self.verbose: print(errstring.format(
                                        checkstring.encode('utf-8'),
                                        self.words[i][j],
                                        self.idxs[i][j]
                                        ))
                                errors['words'][self.idxs[i][j]] = self.words[i][j]
                            tokens = False

                        else:
                            if self.verbose: print(errstring.format(
                                "NO TOKENS",
                                self.words[i][j],
                                self.idxs[i][j]
                                ))
                            errors['words'][self.idxs[i][j]] = self.words[i][j]
        
        # make errors an attribute, if there are any
        self.errors = errors

        # turn self.tokens into LingpyArray
        self.tokens = LingpyArray(self.tokens)

        # turn self.sonars into LingpyArray
        self.alignments = LingpyArray(self.alignments)

        # turn self.sonars into LingPyArray
        self.sonars = LingpyArray(self.sonars)

        # get height
        self.height = len(self.idx)

        # create a dictionary in which data created during the calculation can
        # by stored
        self.data = {}

    def _init_lxs(
            self,
            infile,
            ignore_errors = True
            ):
        """
        Basic routine for input files in lxs-format.
        """
        # load the data
        txt = array(loadtxt(infile),dtype="str")

        # check whether there are cognates in the data
        if 'COG' in txt[0]:
            cog_idx = [i for i in range(txt.shape[1]) if txt[0][i] == 'COG']
            etr_idx = [i-1 for i in cog_idx]
            self.cogs = array(txt[1:,cog_idx],dtype='int')
            self.words = txt[1:,etr_idx]
        else:
            etr_idx = [i for i in range(txt.shape[1]) if txt[0][i].lower() not in
                    ['number','number','words','id']]
            self.words = txt[1:,etr_idx]
            self.cogs = zeros(self.words.shape,dtype='int')
        
        # store the numbers
        self.numbers = array(txt[1:,0],dtype='int')
        
        # check whether headwords are in the data
        if txt[0][1].lower() == 'words':
            self.items = txt[1:,1]
        else:
            self.items = array(self.numbers,dtype='str')

        # create the synonym-idx
        self.idx = {}
        for i,etr in enumerate(self.words):
            try:
                self.idx[self.numbers[i]] += [i]
            except:
                self.idx[self.numbers[i]] = [i]

        # get the language names
        self.taxa = txt[0,etr_idx]
        
        # create a specific format string in order to receive taxa of equal
        # length
        mtax = max([len(t) for t in self.taxa])
        self._txf = '{0:.<'+str(mtax)+'}'

        # get the prosodic strings and the context information
        self.tokens = self.words.copy().tolist()
        self.prostrings = self.words.copy()
        self.sonars = self.words.copy().tolist()
        for i in range(self.words.shape[0]):
            for j in range(self.words.shape[1]):
                if self.words[i][j] != '-':
                    try:
                        tokens = ipa2tokens(self.words[i][j])
                        sonar = [int(x) for x in tokens2class(tokens,art)]
                        prostring = prosodic_string(sonar)
                        self.tokens[i][j] = tokens
                        self.prostrings[i][j] = prostring
                        self.sonars[i][j] = sonar
                    except:
                        if self.verbose: print("[i] Some error with string <"+self.words[i][j]+">")
        
        # turn self.tokens into LingpyArray
        self.tokens = LingpyArray(self.tokens)

        # create an index for all sequences
        self.idxs = self.words.copy()
        self.idxd = {}
        count = 1
        for i,line in enumerate(self.words):
            for j,entry in enumerate(line):
                if entry != '-':
                    self.idxs[i][j] = count
                    self.idxd[count] = (i,j)
                    count += 1
                else:
                    self.idxs[i][j] = 0
        self.idxs = array(self.idxs,dtype='int')
                    
        # get height and width
        self.width = len(self.taxa)
        self.height = len(self.idx)

        # create a dictionary in which data created during the calculation can
        # by stored
        self.data = {}

        # create the self._ifs file for consistency with new csv input format
        # operations (also for output, etc.)
        self._ifs = {
                "ids":(0,'ID'),
                "taxa":(1,'Taxon'),
                "items":(2,'Gloss'),
                "numbers":(3,"GlossID"),
                "words":(4,"IPA"),
                "tokens":(5,"Tokens"),
                }

        # create the self.dataById dictionary for consistency with new csv
        # input format
        self.dataById = {}
        for i in range(self.idxs.shape[0]):
            for j in range(self.idxs.shape[1]):
                if self.idxs[i][j] != 0:
                    self.dataById[self.idxs[i][j]] = [
                            str(self.idxs[i][j]),
                            self.taxa[j],
                            self.items[i],
                            str(self.numbers[i]),
                            self.words[i][j],
                            ' '.join(self.tokens[i][j]).encode('utf-8')
                            ]
        self.dataById[0] = ["ID","Taxon","Gloss","GlossID","IPA","Tokens"]

    def _flatten(
            self,
            idx
            ):
        """
        Return a flat representation of the data in a given semantic slot.
        """

        return [i for i in self.idxs[self.idx[idx]].flatten() if i > 0]
    
    def __getitem__(
            self,
            idx
            ):
        """
        Function returns a specified value for a general ID of the entry.
        """
        
        try:
            data = idx[1]
            idx = abs(idx[0])
        except:
            data = 'w'
            idx = abs(idx)
        
        if data == 'w':
            return self.words[self.idxd[idx]]
        elif data == 'n':
            return self._nbrs[self.idxd[idx]]
        elif data == 'N':
            return self.numbers[self.idxd[idx][0]]
        elif data == 'W':
            return self.weights[self.idxd[idx]]
        elif data == 'r':
            return self.restrictions[self.idxd[idx]]
        elif data == 'c':
            return self.classes[self.idxd[idx]]
        elif data == 'C':
            return self.cogs[self.idxd[idx]]
        elif data == 'l':
            return self.taxa[self.idxd[idx][1]]
        elif data == 'p':
            return self.prostrings[self.idxd[idx]]
        elif data == 'i':
            return self.items[self.idxd[idx][0]]
        elif data == 't':
            return self.tokens[self.idxd[idx]]
        elif data == 'a':
            return self.alignments[self.idxd[idx]]

    def _get_pairs(
            self,
            idxA,
            idxB,
            data,
            repeats = False
            ):
        """
        Return a list of tuples consisting of all pairwise items of the given
        datatype.
        """

        if not repeats:
            try:
                out = self.data[str(idxA)+'/'+str(idxB)+'/'+data+'/f']
            except:
                out = [(self[a,data],self[b,data]) for (a,b) in
                        self._pairs[idxA,idxB] if a > 0]
                self.data[str(idxA)+'/'+str(idxB)+'/'+data+'/f'] = out
        else:
            try:
                out = self.data[str(idxA)+'/'+str(idxB)+'/'+data+'/t']
            except:
                out = [(self[a,data],self[b,data]) for (a,b) in 
                        self._pairs[idxA,idxB] if a != 0]
                self.data[str(idxA)+'/'+str(idxB)+'/'+data+'/t'] = out

        return out

    def __len__(self):
        """
        Return the length of the dataset in terms of the total number of words.
        """
        
        return max(self.idxd)

    def _renumber(
            self,
            loans=False
            ):
        """
        Make the numbers in the data regular.
        @todo: check whether this really works!
        """
        
        num = 1

        for key in self.idx.keys():
            flats = self._flatten(key)
            cogs = [abs(self.cogs[self.idxd[f]]) for f in flats]
            
            uniques = list(set(cogs))
            
            tmp = dict(zip(uniques,range(num,num+len(uniques)+1)))

            for f in flats:
                c = self.cogs[self.idxd[f]]
                if c > 0:
                    self.cogs[self.idxd[f]] = tmp[c]
                elif c < 0:
                    if loans:
                        self.cogs[self.idxd[f]] = -tmp[abs(c)]
                    else:
                        self.cogs[self.idxd[f]] = tmp[abs(c)]

            num += len(uniques)

    def _set_model(
            self,
            model = 'sca',
            merge_vowels = True
            ):
        """
        Define a sequence model for the current calculation and calculate
        several statistics, such as the letter frequency.
        """
        
        # store model and scoring dictionary as attributes
        try:
            self.model = eval(model)
        except:
            self.model = model

        # define the arrays for sound classes and numbers
        self.classes = self.words.copy()
        self._nbrs = self.words.tolist()
        
        # iterate over all entries and change them according to the model
        for i,line in enumerate(self.classes):
            for j,cls in enumerate(line):
                if cls != '-':
                    # check for possible errors
                    try:
                        classes = tokens2class(self.tokens[i][j],self.model)
                        self.classes[i][j] = classes
                        self._nbrs[i][j] = [
                                str(j) + '.' + '.'.join(k) for k in zip(
                                    self.prostrings[i][j],
                                    classes
                                    )
                                    ]
                    except:
                        self.errors['words'][self.idxs[i][j]] = self.words[i][j]

        # convert clss and nbrs into lingpy-array objects
        self._nbrs = LingpyArray(self._nbrs)

        # create the dictionary which stores the frequencies
        self.freqs = {}

        # iterate over all sequences
        for i in range(self.width):
            self.freqs[i] = {}
            for nbr in self._nbrs[:,i]:
                if nbr != '-':
                    for n in nbr:
                        try:
                            self.freqs[i][n] += 1.0
                        except KeyError:
                            self.freqs[i][n] = 1.0

        # create an empty scorer
        self.scorer = {}
        self.score_dict = {}
        
        # iterate over all chars in self._nbrs in order to get all possible
        # combinations of characters and define them as dictionary pairs
        for i in range(self.width):
            for j in range(self.width):
                if i <= j:
                    for charA in list(self.freqs[i].keys())+[str(i)+'.X.-']:
                        for charB in list(self.freqs[j].keys())+[str(j)+'.X.-']:
                            self.scorer[charA,charB] = 0.0
                            try:
                                self.score_dict[charA,charB] = \
                                        self.model.scorer[
                                                charA.split('.')[2],
                                                charB.split('.')[2]
                                                ]
                            except:
                                pass

    def _set_val(
            self,
            scale = None, #(1.2,1.0,1.1),
            factor = 0.3,
            #gop = -3,
            #gep_scale = 0.6,
            restricted_chars = 'T_',
            pairwise_threshold = 0.7
            ):
        """
        Determine the settings for the calculation.
        """

        self.scale = scale
        self.factor = factor
        #self.gop = gop
        #self.gep_scale = gep_scale
        self.restricted_chars = restricted_chars
        self.pairwise_threshold = pairwise_threshold

        # create the weights and the restrictions
        self.restrictions = self.classes.tolist()
        self.weights = self.classes.tolist()

        # iterate over all entries and change them according to the model
        for i,line in enumerate(self.prostrings):
            for j,prostring in enumerate(line):
                if prostring != '-':
                    res = []
                    for k,char in enumerate(prostring):
                        if char in self.restricted_chars:
                            res.append(-(k+1))
                        else:
                            res.append(k+1)
                    self.restrictions[i][j] = res
                    weights = array(
                            prosodic_weights(
                                prostring,
                                scale,
                                factor
                                )
                            )
                        
                    self.weights[i][j] = weights
        
        self.restrictions = LingpyArray(self.restrictions)
        self.weights = LingpyArray(self.weights)

    def _make_pair(
            self,
            idxA,
            idxB
            ):
        """
        Get list of language pairs which serve as input for the align
        algorithm.

        @todo: think of a solution of excluding stuff for the distribution
        calculation and using it later anyway!
        """
        # create the list entry for the dictionary
        self._pairs[idxA,idxB] = []

        # carry out a preprocessing of the word lists in order to avoid that
        # identical words show up in the data
        clsA = self.words[:,idxA].copy()
        clsB = self.words[:,idxB].copy()

        # if the lengths of lists and sets are different, there are duplicate
        # characters which should be eliminated. This is a very rough procedure
        # which may well lead to a loss of information
        if len(clsA) != len(set(clsA)):
            tmp = []
            for i,cls in enumerate(clsA):
                if cls in tmp and cls != '-':
                    clsA[i] = '--'
                else:
                    tmp.append(cls)
        if len(clsB) != len(set(clsB)):
            tmp = []
            for i,cls in enumerate(clsB):
                if cls in tmp and cls != '-':
                    clsB[i] = '--'
                else:
                    tmp.append(cls)

        # fill the lists
        for key in self.idx.keys():
            # get all lists for the language pairs
            listA = clsA[self.idx[key]]
            listB = clsB[self.idx[key]]
            for i,lA in enumerate(listA):
                for j,lB in enumerate(listB):
                    if '-' not in (lA,lB):
                        
                        # get the pairs
                        pairA = self.idxs[self.idx[key],idxA][i]
                        pairB = self.idxs[self.idx[key],idxB][j]

                        if '--' not in (lA,lB):
                            self._pairs[idxA,idxB].append((pairA,pairB))
                        else:
                            self._pairs[idxA,idxB].append((-pairA,-pairB))
            
    def _make_all_pairs(
            self
            ):
        """
        Create pairwise lists for all languages.
        """

        self._pairs = {}
        self._alignments = {}
        self.scores = {}

        # iterate over all languages
        for i in range(self.width):
            for j in range(self.width):
                if i <= j:
                    self._make_pair(i,j)

    def _scale(
            self,
            x,
            factor = 0.01
            ):
        """
        Scaling factor for distances in dependence of sequence length.
        """
        try:
            return self.scales[x-1]
        except:
            s = 0
            self.scales = [1]
            for i in range(1,x):
                s += (x / i) * (factor ** i) * (1 - factor) ** (x-i)
                self.scales.append(1+s)
            return self.scales[x-1]

    def _self_score(
            self,
            seq,
            scorer
            ):
        """
        Return the score for an alignment with itself.
        """
        score = 0.0
        new_seq = [i for i in seq if '-' not in i]
        
        for a in new_seq:
            score += scorer[a,a]

        return score

    def _pair_score(
            self,
            seqA,
            seqB,
            scorer
            ):
        """
        Return the score for an alignment of two sequences.
        """
        score = 0.0
        for a,b in zip(seqA,seqB):
            if '-' not in (a,b):
                score += scorer[a,b]
            else:
                score += -1.0 

        return score

    def _distance_score(
            self,
            almA,
            almB,
            scoreAB,
            scorer
            ):
        """
        Function calculates the Downey et al. (2008) distance score for
        pairwise sequence alignments.
        """
        
        # calculate the self_scores
        scoreA = self._self_score(almA,scorer)
        scoreB = self._self_score(almB,scorer)
        
        try:
            # @check@ whether it changes the results if the scores are normalized
            # to not exceed 1.0
            score = (1.0 - (2.0 * scoreAB / (scoreA + scoreB))) 
            if score > 1.0:
                return score
            else:
                return score

        except ZeroDivisionError:   
            #if self.verbose: print(almA,almB,scoreAB)
            if self.verbose: print("[!] Warning, float division!")
            return 10.0
   
    def _align_pairwise(
            self,
            idxA,
            idxB,
            mode = ('global',-2,0.5),
            ):
        """
        Align all words of two languages pairwise.
        """
                
        if idxA == idxB:
            numbers = self._get_pairs(idxA,idxA,'n')
            alignments = [
                    (
                        a[0],
                        a[0],
                        self._self_score(a[0],self.score_dict)
                        )
                    for a in numbers]
        else:
            numbers = self._get_pairs(idxA,idxB,'n')
            prostrings = self._get_pairs(idxA,idxB,'p')
            restrictions = self._get_pairs(idxA,idxB,'r')
            weights = [(list(mode[1] * a),list(mode[1] * b)) for a,b in 
                self._get_pairs(idxA,idxB,'W')]

            alignments = align_sequence_pairs(
                    numbers,
                    weights,
                    restrictions,
                    prostrings,
                    self.score_dict,
                    mode[2],
                    self.factor,
                    mode[0]
                    )
        
        # change alms if mode is local
        # @todo: in order to avoid loosing information, this part should
        # eventually be changed in such a way that all the rest of sequence
        # parts in local alignments is reflected as gaps, otherwise, we can
        # simply think of using dialign instead of local analyses
        if mode[0] == 'local':
            for i,alm in enumerate(alignments):
                almA = [k for k in alm[0] if k != '*']
                almB = [k for k in alm[1] if k != '*']
                score = alm[2]
                alignments[i] = (almA,almB,score)
        
        return alignments
    
    def _get_correspondences(
            self,
            alignments,
            idxA,
            idxB,
            heuristic = 'distance'
            ):
        """
        Function returns correspondence statistics. The threshold determines
        the maximum distance between two sequences which is included in the
        overall score.
        """
        reg_dist = {}

        if idxA != idxB and heuristic == 'distance':
            for almA,almB,scoreAB in alignments:
                
                score = self._distance_score(almA,almB,scoreAB,self.score_dict)
                if score <= self.pairwise_threshold:
                    for a,b in zip(almA,almB):
                        try:
                            reg_dist[a,b] += 1.0
                        except KeyError:
                            reg_dist[a,b] = 1.0
        
        # a new attempt, we use a preprocessing and cognate clustering with
        # simple sca and then only consider those words found by the rough
        # algorithm for the calculation of sound correspondences
        elif heuristic == 'sca':
            for i,(almA,almB,scoreAB) in enumerate(alignments):
                x,y = self._pairs[idxA,idxB][i]
                # remember that we excluded pairs from the calc which occur
                # more than ones, i.e. identical strings, so we follow the same
                # procedure here
                if x > 0 and y > 0:
                    x = self.cogs[self.idxd[x]]
                    y = self.cogs[self.idxd[y]]
                    if x == y:
                        for a,b in zip(almA,almB):
                            try:
                                reg_dist[a,b] += 1.0
                            except:
                                reg_dist[a,b] = 1.0

        elif idxA == idxB:
            counter = 0
            for almA,almB,scoreAB in alignments:
                
                score = self._distance_score(almA,almB,scoreAB,self.score_dict)

                if score == 0:
                    for a,b in zip(almA,almB):
                        try:
                            reg_dist[a,b] += 1.0
                        except KeyError:
                            reg_dist[a,b] = 1.0
                    counter += 1

        # change gap notation
        for a,b in list(reg_dist.keys()):
            
            # change for gap positions
            if a == '-':
                reg_dist[str(idxA)+'.X.-',b] = reg_dist[a,b]
            elif b == '-':
                reg_dist[a,str(idxB)+'.X.-'] = reg_dist[a,b]

        return reg_dist

    def _random_align_pairwise(
            self,
            idxA,
            idxB,
            mode = ('global',-2,0.5),
            runs = 50
            ):
        """
        Align sequences pairwise, whereas all lists are shuffled in order to
        get a random distribution of the data.
        """
        numbers = self._get_pairs(idxA,idxB,'n')
        prostrings = self._get_pairs(idxA,idxB,'p')
        restrictions = self._get_pairs(idxA,idxB,'r')
        weights = [(list(mode[1] * a),list(mode[1] * b)) for a,b in 
                self._get_pairs(idxA,idxB,'W')]
       
        rand_dist = random_align_sequence_pairs(
                numbers,
                weights,
                restrictions,
                prostrings,
                self.score_dict,
                mode[2],
                self.factor,
                mode[0],
                runs
                )
        
        # change gap symbols in order to make them comparable to pairwise
        # notation
        for a,b in list(rand_dist.keys()):
            if a == '-':
                rand_dist[str(idxA)+'.X.-',b] = rand_dist[a,b]
            elif b == '-':
                rand_dist[a,str(idxB)+'.X.-'] = rand_dist[a,b]

        return rand_dist

    def _join_dist(
            self,
            dists
            ):
        """
        Function joins two or more distributions by averaging them.
        """
        
        if len(dists) == 1:
            return dists[0]
        
        out_dist = {}
        
        keys = []
        for dist in dists:
            keys += list(dist.keys())
        
        keys = set(keys)

        for key in keys:
            vals = []
            for dist in dists:
                try:
                    vals.append(dist[key])
                except:
                    vals.append(0.0)

            out_dist[key] = sum(vals) / len(dists)

        return out_dist

    def _expand_scorer_pairwise(
            self,
            idxA,
            idxB,
            runs = 50,
            modes = (('global',-2,0.5),('local',-1,0.5)),
            ratio = (3,2),
            heuristic = 'distance'
            ):
        """
        Add new scores to the global scoring dictionary for the dataset. Scores
        are determined by taking the logarithm of the division of the square of
        attested and expected frequencies for each character combination.
        """
        
        # get distribution for random alignments
        expected = []
        for calc in modes:
            exp = self._random_align_pairwise(
                    idxA,
                    idxB,
                    calc,
                    runs
                    )
            expected.append(exp)
        expected = self._join_dist(expected)
        
        # get distribution for real alignments
        attested = []
        for calc in modes:
            att = self._get_correspondences(
                    self._align_pairwise(
                        idxA,
                        idxB,
                        calc,
                        ),
                    idxA,
                    idxB,
                    heuristic
                    )
            attested.append(att)
        attested = self._join_dist(attested)
        
        # determine an average gop for the cal
        gop = sum([m[2] for m in modes]) / len(modes)

        # update the scorer
        for charA in list(self.freqs[idxA].keys())+[str(idxA)+'.X.-']:
            for charB in list(self.freqs[idxB].keys())+[str(idxB)+'.X.-']:
                try:
                    exp = expected[charA,charB]
                except KeyError:
                    exp = False
                try:
                    att = attested[charA,charB]
                except KeyError:
                    att = False
                
                # if the residue pair is only attested once in the dataset,
                # this may well be a coincidence. Therefore, the score is
                # set to 0.01, in order to avoid that possible coincidences bear to
                # much weight. note that this value eventually has to be
                # modified and further investigated.
                # XXX note especially that this should not be done for
                # identical language pairs being compared, since here,
                # everything should be indicative XXX
                if att <= 1 and idxA != idxB:
                    att = False
                
                # if there are values for both attested and expected residue
                # pairs, the algorithm follows Kessler (2000) in so far, as the
                # values are squared in order to make them more "indicative",
                # furthermore, the binary logarithm of the division of the
                # attested and the expected values is taken in order to
                # retrieve a dictionary score
                if att and exp:
                    score = log2((att ** 2) / (exp ** 2))

                # if a residue pair is only attested and not expected (which is
                # possible, if the number of shuffled iterations is too low),
                # we simply assume that the square of the expected value is
                # 0.01. this might result in a certain bias which should be
                # kept in mind when using the algorithm
                elif att and not exp:
                    score = log2((att ** 2) / 0.01)

                # if a residue pair is only expected but not attested, this
                # certainly should result in a negative values. in order to
                # avoid problematic calculations, we simply set the value to
                # -5. this may, again result in a certain bias which should be
                # kept in mind.
                elif exp and not att:
                    score = -5.0 

                # if a residue pair is neither expected nor attested, the score
                # should surely be very low, and we simply set it to -90 in
                # order to avoid the algorithm to match such residues during
                # the alignment process.
                elif not exp and not att:
                    score = -90.0
                
                # get the scores for the regular scoring dictionary. these
                # scores are combined with the correspondence-based scoring
                # scheme in order to cope for the fact that information
                # available in small word lists might be low. the combination
                # is based on a certain ratio by which both values are
                # combined. the current default is 2 : 1, i.e. the
                # correspondence-based scoring scheme counts twice as much as
                # the regular scoring scheme of the alignment algorithm being
                # applied. in case of a gap, the regular gap score as it is
                # defined in the beginning is used to represent the regular
                # scoring function.
                if '-' not in charA+charB:
                    sim = self.score_dict[charA,charB]
                else:
                    sim = gop
                
                # combine the scores according to the given ratio
                # give only half of the scores to vowels, since they otherwise
                # weigh too heavy
                this_score = (ratio[0] * score + ratio[1] * sim) / sum(ratio)

                if charA[2] in 'Vv<>' and charB[2] in 'Vv<>':
                    self.scorer[charA,charB] = self.vscale * this_score
                elif charA[2] in 'T_' and charB[2] in 'T_':
                    self.scorer[charA,charB] = self.vscale * this_score
                # XXX add a negative starter for gaps
                #elif 'X' in charA[2]+charB[2]:
                #    self.scorer[charA,charB] = this_score - 1
                else:
                    self.scorer[charA,charB] = this_score
                #self.scorer[
                #        charA,
                #        charB
                #        ] = (ratio[0] * score + ratio[1] * sim) / sum(ratio)

    def _expand_scorer(
            self,
            runs = 50,
            modes = (('global',-3,0.5),('local',-1,0.5)),
            ratio = (1,1),
            heuristic = 'distance'
            ):
        """
        Carry out pairwise and random alignments and calculate new scores for
        the library scorer.
        """
        for i in range(self.width):
            for j in range(self.width):
                if i <= j:
                    if self.verbose: print("[i] Calculating scores for",self.taxa[i],"and",self.taxa[j],"...")
                    self._expand_scorer_pairwise(i,j,runs,modes,ratio,heuristic)

    def _get_pairwise_scores(
            self,
            idxA,
            idxB,
            mode = 'overlap',
            scale = None, #(1.0,1.0,1.0),
            factor = 0.3,
            gep_scale = 0.5,
            score_mode = 'library',
            gop = -2,
            ):
        """
        Function calculates distance scores for pairwise alignments of the
        wordlists of two languages.
        """

        if score_mode == 'library':
            scorer = self.scorer

            # calculate alignments
            # determine weights on the basis of most probable gaps
            weights = []
            numbers = self._get_pairs(idxA,idxB,'n',True)
            for a,b in numbers:
                wA,wB = [],[]
                for n in a:
                    wA.append(self.scorer[n,str(idxB)+'.X.-'])
                for n in b:
                    wB.append(self.scorer[str(idxA)+'.X.-',n])
                weights.append((wA,wB))
            restrictions = self._get_pairs(idxA,idxB,'r',True)
            prostrings = self._get_pairs(idxA,idxB,'p',True)

            # carry out alignments
            alignments = align_sequence_pairs(
                    numbers,
                    weights,
                    restrictions,
                    prostrings,
                    self.scorer,
                    gep_scale,
                    factor,
                    mode
                    )

        # simple score_mode
        elif score_mode == 'sca':
            scorer = self.score_dict
            numbers = self._get_pairs(idxA,idxB,'n',True)
            restrictions = self._get_pairs(idxA,idxB,'r',True)
            prostrings = self._get_pairs(idxA,idxB,'p',True)
            weights = [
                    (
                        list(gop * array(prosodic_weights(
                            a,
                            scale,
                            factor
                            ))),
                        list(gop * array(prosodic_weights(
                            b,
                            scale,
                            factor
                            )))
                        ) for (a,b) in prostrings]

            alignments = align_sequence_pairs(
                    numbers,
                    weights,
                    restrictions,
                    prostrings,
                    self.score_dict,
                    gep_scale,
                    factor,
                    mode
                    )

        # turchin score-mode
        elif score_mode == 'turchin':
            for i,(a,b) in enumerate(self._pairs[idxA,idxB]):
                
                a,b = abs(a),abs(b)
                tmpA = tokens2class(
                        self.tokens[self.idxd[a]],
                        dolgo
                        ).replace('V','')
                tmpB = tokens2class(
                        self.tokens[self.idxd[b]],
                        dolgo
                        ).replace('V','')
                
                if tmpA[0:2] == tmpB[0:2]:
                    dist = 0.0
                else:
                    dist = 1.0

                self.scores[a,b] = dist

        # edit distance
        elif score_mode == 'edit-dist':
            for i,(a,b) in enumerate(self._pairs[idxA,idxB]):

                a,b = abs(a),abs(b)
                tmpA = list(unicode(self[a],'utf-8'))
                tmpB = list(unicode(self[b],'utf-8'))

                dist = edit_dist(tmpA,tmpB)

                self.scores[a,b] = dist
        
        elif score_mode == 'edit-tokens':
            for i,(a,b) in enumerate(self._pairs[idxA,idxB]):

                a,b = abs(a),abs(b)
                tmpA = self[a,'t']
                tmpB = self[b,'t']

                dist = edit_dist(tmpA,tmpB)

                self.scores[a,b] = dist

        # change alms if mode is local
        if mode == 'local':
            for i,alm in enumerate(alignments):
                almA = [k for k in alm[0] if k != '*']
                almB = [k for k in alm[1] if k != '*']
                score = alm[2]
                alignments[i] = (almA,almB,score)

        if score_mode in ['sca','library']:
            for i,(almA,almB,scoreAB) in enumerate(alignments):
                
                # get the pairs
                pairA,pairB = self._pairs[idxA,idxB][i]
                pairA,pairB = abs(pairA),abs(pairB)
                
                # store the alignments
                self._alignments[pairA,pairB] = [almA,almB,scoreAB]
                
                # calculate the distance
                distAB = self._distance_score(almA,almB,scoreAB,scorer)

                self.scores[pairA,pairB] = distAB

    def _get_all_pairwise_scores(
            self,
            mode = 'overlap',
            scale = None, #(1.0,1.0,1.0),
            factor = 0.3,
            gep_scale = 0.5,
            score_mode = 'library',
            gop = -2,
            ):
        """
        Calculate all pairwise scores for the current dataset.
        """

        for i in range(self.width):
            for j in range(self.width):
                if i <= j:
                    self._get_pairwise_scores(
                            i,
                            j,
                            mode,
                            scale,
                            factor,
                            gep_scale,
                            score_mode,
                            gop
                                )
    def _align_profile(
            self,
            almsA,
            almsB,
            sonA,
            sonB,
            return_self = False
            ):
        """
        Function is a simplified MSA-version similar, but - due to the task -
        not identical to proper SCA alignment.
        """
        # turn the alignments into profiles
        profileA = array(almsA).transpose()
        profileB = array(almsB).transpose()

        # get the prosodic profiles
        prosA = prosodic_string(sonA)
        prosB = prosodic_string(sonB)
        
        # get the languages involved
        tmpA,tmpB = [],[]
        for x in profileA:
            tmpA += [y[0] for y in x if y != 'X']
        for x in profileB:
            tmpB += [y[0] for y in x if y != 'X']

        numsA = list(set(tmpA))
        numsB = list(set(tmpB))

        m,o = profileA.shape
        n,p = profileB.shape

        lstA = list(range(1,m+1))
        lstB = list(range(1,n+1))

        # get the restricted chars
        resA,resB = [l for l in lstA],[l for l in lstB]
        for i,a in enumerate(prosA):
            if a in self.restricted_chars:
                resA[i] = -(lstA[i])
        for i,a in enumerate(prosB):
            if a in self.restricted_chars:
                resB[i] = -(lstB[i])
        
        # calculate the weights for all positions by summing up the individual
        # scores for the segment pairs. currently, the scores are slightly
        # lowered, since they would otherwise be too high and lower the recall
        weightsA,weightsB = [],[]
        
        for col in profileA:
            w = 0.0
            c = 0.0
            for cell in [k for k in col if k != 'X']:
                for x in numsA:
                    try:
                        w += self.scorer[cell,x+'.X.-']
                        c += 1.0
                    except:
                        w += self.scorer[x+'.X.-',cell]
                        c += 1.0
            weightsA.append((w / c) * 0.5 ) 

        for col in profileB:
            w = 0.0
            c = 0.0
            for cell in [k for k in col if k != 'X']:
                for x in numsB:
                    try:
                        w += self.scorer[cell,x+'.X.-']
                        c += 1.0
                    except:
                        w += self.scorer[x+'.X.-',cell]
                        c += 1.0
            weightsB.append((w / c) * 0.5) 

        weightsA = array(weightsA) * prosodic_weights(prosA)
        weightsB = array(weightsB) * prosodic_weights(prosB)
        weightsA = weightsA.tolist()
        weightsB = weightsB.tolist()

        # create the scoring dictionary
        tmp = {}
        for I in lstA:
            for J in lstB:
                score = 0.0
                counter = 0.0
                if prosA[I-1] != prosB[J-1]:
                    match = False
                else:
                    match = True
                for charA in profileA[I-1]:
                    for charB in profileB[J-1]:
                        if 'X' not in (charA,charB):
                            try:
                                this_score = self.scorer[charA,charB]
                            except KeyError:
                                this_score = self.scorer[charB,charA]
                            if match:
                                this_score += this_score * self.factor
                            score += this_score #self.scorer[charA,charB]
                            counter += 1.0

                        else:
                            counter += 0.0
                tmp[I,J] = score / counter

        alignedA,alignedB,sim = align_pairwise(
                lstA,
                lstB,
                weightsA,
                weightsB,
                resA,
                resB,
                len(lstA) * 'X',
                len(lstB) * 'X',
                tmp,
                0.5,
                0.0,
                'overlap'
                )

        if return_self:
            return sim
        
        new_sonA = [s for s in sonA]
        new_sonB = [s for s in sonB]
        profileA = profileA.tolist()
        profileB = profileB.tolist()
        for i in range(len(alignedA)):
            if alignedA[i] == '-':
                profileA.insert(i,o *['X'])
                new_sonA.insert(i,0)
            elif alignedB[i] == '-':
                profileB.insert(i,p * ['X'])
                new_sonB.insert(i,0)

        profileA = array(profileA).transpose().tolist()
        profileB = array(profileB).transpose().tolist()

        # get the consensus of the sonority profile
        cons = [int(i[i!=0].mean().round()) for i in
                array([new_sonA,new_sonB]).transpose()]

        return profileA+profileB,cons,sim
        
    def _alm_cluster(
            self,
            idx,
            threshold=0.5
            ):
        """
        Flate clustering method which is based on simultaneous alignment and
        scoring of sequences along a guide tree until a threshold of distance
        scores is reached.
        """
        # get the flats
        flats = self._flatten(idx)

        # create the cluster dictionary
        clusters = dict([(i,[i]) for i in range(len(flats))])

        # create a dictionary for the alignments
        alm_dict = dict([(i,[self._nbrs[self.idxd[flats[i]]]]) for i in
            range(len(flats))])

        # create a dictionary for the sonority-information
        son_dict = dict([(i,self.sonars[self.idxd[flats[i]]]) for i in
            range(len(flats))])
        
        # create the matrix
        matrix = []
        tree = []

        # fill in the matrix with the calculated distance scores
        for i,idxA in enumerate(flats):
            for j,idxB in enumerate(flats):
                if i < j:
                    try:
                        matrix.append(self.scores[idxA,idxB])
                    except:
                        matrix.append(self.scores[idxB,idxA])

        # turn the flat matrix into a redundant matrix
        matrix = squareform(matrix)

        # cluster the data
        _upgma(clusters,matrix,tree)
        
        stop = len(tree) - 1
        newick = dict([(i,[i]) for i in range(len(flats))])
        x = len(flats)
        clusters = {}
        alms = {}
        y = 1
        newick[x] = True
        # align the seqs along the guide tree until a limit is reached
        for i,(a,b,c,d) in enumerate(tree):
            if newick[a] and newick[b]:
                alm,son,sAB = self._align_profile(
                        alm_dict[a],
                        alm_dict[b],
                        son_dict[a],
                        son_dict[b]
                        )
                sA = self._align_profile(
                        alm_dict[a],
                        alm_dict[a],
                        son_dict[a],
                        son_dict[a],
                        True
                        )
                sB = self._align_profile(
                        alm_dict[b],
                        alm_dict[b],
                        son_dict[b],
                        son_dict[b],
                        True
                        )
                try:
                    score = 1 - (2 * sAB) / (sA + sB)
                except ZeroDivisionError:
                    score = 10

                if score <= threshold and i != stop:
                    alm_dict[x+i] = alm
                    newick[x+i] = newick[a]+newick[b]
                    son_dict[x+i] = son
                elif score > threshold and i != stop:
                    for A,B in zip(newick[a],alm_dict[a]):
                        clusters[A] = y
                        alms[A] = B
                    y += 1
                    for A,B in zip(newick[b],alm_dict[b]):
                        clusters[A] = y
                        alms[A] = B
                    y += 1
                    newick[x+i] = False
                elif score <= threshold and i == stop:
                    for A,B in zip(newick[a] + newick[b],alm):
                        clusters[A] = y
                        alms[A] = B
                    newick[x+i] = newick[a]+newick[b]
                elif score > threshold and i == stop:
                    for A,B in zip(newick[a],alm_dict[a]):
                        clusters[A] = y
                        alms[A] = B
                    y += 1
                    for A,B in zip(newick[b],alm_dict[b]):
                        clusters[A] = y
                        alms[A] = B
            else:
                newick[x+i] = False
                if newick[a]:
                    for A,B in zip(newick[a],alm_dict[a]):
                        clusters[A] = y
                        alms[A] = B
                    y += 1
                if newick[b]:
                    for A,B in zip(newick[b],alm_dict[b]):
                        clusters[A] = y
                        alms[A] = B
                    y += 1
        
        if not clusters:
            out = dict([(i,[i]) for i in range(x)])
        else:
            out = {}
            for i,j in clusters.items():
                try:
                    out[j] += [i]
                except KeyError:
                    out[j] = [i]
        # get the keys for the clusters
        count = 1
        for key in out:
            for val in out[key]:
                flats[val] = (count,flats[val])
            count += 1
        
        # store the alignments in tokenized form
        for key in alms.keys():
            a,b = flats[key]
            for i,v in enumerate(alms[key]):
                if v == 'X':
                    self.alignments[self.idxd[b]].insert(i,'-') 
        
        return flats

    def _cluster(
            self,
            idx,
            threshold = 0.5,
            debug = True
            ):
        """
        Cluster the data in order to get unified cognate judgments throughout
        the word lists.
        """
        
        # get the flats
        flats = self._flatten(idx)
        
        # create cluster dictionary
        clusters = dict([(i,[i]) for i in range(len(flats))])


        # create the matrix
        matrix = []
        
        # fill in the matrix with the calculated distance scores
        for i,idxA in enumerate(flats):
            for j,idxB in enumerate(flats):
                if i < j:
                    try:
                        matrix.append(self.scores[idxA,idxB])
                    except:
                        matrix.append(self.scores[idxB,idxA])
        
        # turn the flat matrix into a redundant matrix
        matrix = squareform(matrix)

        if debug:
            for line in matrix:
                print(" ".join(["{0:2f}".format(l) for l in line]))

        # cluster the data
        _flat_upgma(clusters,matrix,threshold)

        # flat neighbor XXX test
        #clusters = _flat_neighbor(matrix,range(len(flats)),threshold)

        # get the keys for the clusters
        count = 1
        for key in clusters:
            for val in clusters[key]:
                flats[val] = (count,flats[val])
            count += 1
        
        return flats
    
    def _get_cognates(
            self,
            threshold,
            clr_mode = "simple"
            ):
        """
        Calculate possible cognates from the dataset.
        """
        
        # determine the counter
        count = 0
        
        # iterate over all semantic slots and cluster the data
        if 'cogid' in ','.join(self.dataById[0]).lower():
            for key in self.idx:
                if clr_mode == 'alm':
                    flats = self._alm_cluster(key,threshold)
                else:
                    flats = self._cluster(key,threshold)
                
                for a,b in flats:
                    self.cogs[self.idxd[b]] = a + count
                    self.dataById[self.idxs[self.idxd[b]]][self._ifs['cogs'][0]] = str(a + count)

                count += max([k[0] for k in flats])

        elif 'cogid' not in ','.join(self.dataById[0]).lower():
            for key in self.idx:
                if clr_mode == 'alm':
                    flats = self._alm_cluster(key,threshold)
                else:
                    flats = self._cluster(key,threshold)
                
                for a,b in flats:
                    self.cogs[self.idxd[b]] = a + count
                    
                count += max([k[0] for k in flats])

    def _etym_dict(
            self,
            include_loans = False
            ):
        
        self.etym_dict = {}
        """
        A dictionary which contains the information regarding cognacy of a
        specified dataset in 'etymological' format, i.e. every cognate is given
        a specific ID and the words and meanings corresponding to this ID are
        listed as values.
        """

        # get all ids present in the dataset
        if include_loans:
            dict_ids = list(set([abs(i) for i in self.cogs.flatten() if i != 0]))
        else:
            dict_ids = list(set([i for i in self.cogs.flatten() if i != 0]))

        # iterate over the cognates and append each word corresponding to a
        # given ID to the dictionary
        for dict_id in dict_ids:
            if dict_id != 0:
                self.etym_dict[dict_id] = [
                        ['-' for i in range(len(self.taxa))],
                        ['-' for i in range(len(self.taxa))],
                        [[] for i in range(len(self.taxa))],
                        [lng for lng in self.taxa],
                        [],
                        []
                        ]
        
        for i in range(len(self.cogs)):
            for j in range(len(self.cogs[i])):
                if self.cogs[i][j] != 0:
                    if include_loans:
                        tmp = abs(self.cogs[i][j])
                    else:
                        tmp = self.cogs[i][j]

                    if self.etym_dict[tmp][0][j] != '-':
                        self.etym_dict[tmp][0][j] += \
                            [self.words[i][j]]
                        self.etym_dict[tmp][1][j] += \
                            [self.items[i]]
                        self.etym_dict[tmp][2][j] += \
                            [self.numbers[i]]
                        self.etym_dict[tmp][4].append(self.idxs[i][j])
                        
                        # trace borrowings
                        if self.cogs[i][j] < 0:
                            self.etym_dict[tmp][5].append(1)
                        else:
                            self.etym_dict[tmp][5].append(0)

                    else:
                        self.etym_dict[tmp][0][j] = \
                            [self.words[i][j]]
                        self.etym_dict[tmp][1][j] = \
                            [self.items[i]]
                        self.etym_dict[tmp][2][j] = \
                            [self.numbers[i]]
                        self.etym_dict[tmp][4].append(self.idxs[i][j])

                        # trace borrowings
                        if self.cogs[i][j] < 0:
                            self.etym_dict[tmp][5].append(1)
                        else:
                            self.etym_dict[tmp][5].append(0)

    def pairwise_distances(self):
        """
        Calculate the lexicostatistical distance between all taxa.

        Returns
        ----------
        dist_matrix : `numpy.array`
            A two-dimensional array containing the scores for the pairwise
            distances between all taxa.

        Examples
        --------
        Load the benchmark file :file:`SLV.lxs` which contains manually
        conducted cognate judgments.

        >>> from lingpy import *
        >>> lex = LexStat(get_file('SLV.lxs'))
        >>> dist = lex.pairwise_distances()
        >>> formstring = '{0[0]:.2f} {0[1]:.2f} {0[2]:.2f} {0[3]:.2f}'
        >>> for line in dist: if self.verbose: print(formstring.format(line))
        0.00 0.15 0.18 0.17
        0.15 0.00 0.20 0.10
        0.18 0.20 0.00 0.20
        0.17 0.10 0.20 0.00

        """

        dist_matrix = []

        for i in range(len(self.taxa)):
            for j in range(len(self.taxa)):
                if i < j:
                    # iterate through both lists and store, whether two entries
                    # are the same or not, ignore gapped entries
                    temp = []
                    langA = self.cogs[:,i]
                    langB = self.cogs[:,j]
                    for key in self.idx.keys():
                        ind = self.idx[key]
                        tmp = []
                        if max(langA[ind]) > 0 and max(langB[ind]) > 0:
                            for num in langA[ind]:
                                if num > 0 and num in langB[ind]:
                                    tmp.append(1)
                                else:
                                    tmp.append(0)
                        if len(tmp) != 0:
                            temp.append(max(tmp))
                    hits = sum(temp)
                    counts = len(temp)
                    try:
                        calc = 1 - float(hits) / counts
                    except:
                        print(self.taxa[i],self.taxa[j])
                    dist_matrix.append(calc)
        dist_matrix = squareform(dist_matrix)

        return dist_matrix

    def _make_cognate_pairs(
            self,
            include_loans = False
            ):
        """
        Function returns a dictionary of all word-pairs (id as key) along with
        an indicator regarding cognacy.
        """

        try:
            self._pairs
        except:
            self._set_model()
            self._set_val()
            self._make_all_pairs()

        # iterate over all word-pairs and determine the cognacy
        cognates = {}
        
        for i in range(self.width):
            for j in range(self.width):
                if i <= j:
                    for pairA,pairB in self._pairs[i,j]:
                        pairA,pairB = abs(pairA),abs(pairB)

                        if pairA != pairB:
                            cogA = self[pairA,'C']
                            cogB = self[pairB,'C']

                            if cogA == cogB and cogA > 0:
                                cognates[pairA,pairB] = 1
                            elif cogA < 0 or cogB < 0:
                                if not include_loans:
                                    cognates[pairA,pairB] = 0
                                else:
                                    if abs(cogA) == abs(cogB):
                                        cognates[pairA,pairB] = 1
                                    else:
                                        cognates[pairA,pairB] = 0
                            else:
                                cognates[pairA,pairB] = 0
        
        return cognates

    def analyze(
            self,
            threshold,
            score_mode = 'library',
            model = 'sca',
            merge_vowels = True,
            scale = None, 
            factor = 0.3,
            restricted_chars = 'T_',
            pairwise_threshold = 0.7,
            runs = 100,
            modes = (('global',-2,0.5),('local',-1,0.5)),
            ratio = (3,2),
            mode = 'overlap',
            clr_mode = 'simple',
            vscale = 0.5,
            heuristic = 0.6,
            verbose = False,
            gep_scale = 0.5,
            gop = -2.0
            ):
        """
        Conduct automatic cognate judgments following the method of :evobib:`List2012b`.

        Parameters
        ----------
        threshold : float
            The threshold which is used for the flat cluster analysis.

        score_mode : { 'library', 'sca', 'turchin', 'edit-dist', 'edit-tokens' }
            Define the `score_mode` on which the calculation of pairwise
            distances is based. Select between:

            * 'library' -- the distance scores are based on the
              language-specific scoring schemes as described in
              :evobib:`List2012b` (this is the default),

            * 'sca' -- the distance scores are based on the
              language-independent SCA distance (see :evobib:`List2012b`),

            * 'turchin' -- the distance scores are based on the approach
              described in :evobib:`Turchin2010`,

            * 'edit-dist"' -- the distance scores are based on the normalized
              edit distance (:evobib:`Levenshtein1966`), and

            * 'edit-tokens' -- the distance scores are based on the normalized
              edit distance, yet the scores are derived from the tokenized
              representation of the sequences and not from their raw,
              untokenized form.

        model : string (default="sca")
            A string indicating the name of the
            :py:class:`~lingpy.data.model.Model` object that shall be used
            for the analysis.
            Currently, three models are supported:
            
            * "dolgo" -- a sound-class model based on :evobib:`Dolgopolsky1986`,

            * "sca" -- an extension of the "dolgo" sound-class model based on
              :evobib:`List2012a`, and
            
            * "asjp" -- an independent sound-class model which is based on the
              sound-class model of :evobib:`Brown2008` and the empirical data of
              :evobib:`Brown2011`.

        merge_vowels : bool (default=True)
            Indicate, whether neighboring vowels should be merged into
            diphtongs, or whether they should be kept separated during the
            analysis.

        gop : int (default=-5)
            The gap opening penalty (gop) on which the analysis shall be based.

        gep_scale : float (default=0.6)
            The factor by which the penalty for the extension of gaps (gap
            extension penalty, GEP) shall be decreased. This approach is
            essentially inspired by the extension of the basic alignment
            algorithm for affine gap penalties by :evobib:`Gotoh1982`.

        scale : tuple or list (default=(3,1,2))
            The scaling factors for the modificaton of gap weights. The first
            value corresponds to sites of ascending sonority, the second value
            to sites of maximum sonority, and the third value corresponds to
            sites of decreasing sonority.

        factor : float (default=0.3)
            The factor by which the initial and the descending position shall
            be modified.
       
        restricted_chars : string (default="T")
            Define which characters of the prosodic string of a sequence
            reflect its secondary structure (cf. :evobib:`List2012a`) and should
            therefore be aligned specifically. This defaults to "T", since this
            is the character that represents tones in the prosodic strings of
            sequences.

        pairwise_threshold : float (default=0.7)
            Only those sequence pairs whose distance is beyond this threshold
            will be considered when determining the distribution of attested
            segment pairs.

        runs : int (default=100)
            Define how many times the perturbation method shall be carried out
            in order to retrieve the expected distribution of segment pairs.

        modes : tuple or list (default = ("global","local"))
            Define the alignment modes of the pairwise analyses which are
            carried out in order to create the language-specific scoring scheme.

        ratio : tuple (default=(1,1))
            Define the ratio by which the traditional scoring scheme and the
            correspondence-based scoring scheme contribute to the actual
            library-based scoring scheme. 

        mode : string (default = "overlap")
            Define the alignment mode which is used in order to calculate
            pairwise distance scores from the language-specific scoring
            schemes.

        """
        # define the factor and restricted chars as an attribute
        self.factor = factor
        self.restricted_chars = restricted_chars
        self.verbose = verbose

        # in order to save time, the values for a given calculation are dumped
        # into a binary file and loaded, if the calculation already has been
        # carried out before
        vals = {
                'model':repr(model),
                'merge_vowels':str(merge_vowels),
                'factor':str(factor),
                'mgop':'gop:'+'-'.join([str(m[1]) for m in modes]),
                'restricted_chars':restricted_chars,
                'GEP':'gep:'+'-'.join(['{0:.2f}'.format(m[2]) for m in
                    modes]),
                'GOP' : str(gop),
                'gep_scale':str(gep_scale),
                'modes':'-'.join([m[0] for m in modes]),
                'ratio':str(ratio[0])+'-'+str(ratio[1]),
                'pairwise_threshold':str(pairwise_threshold),
                'runs':str(runs),
                'heuristic':str(heuristic),
                'vscale' : '{0:.1f}'.format(vscale)
                }

        self.vscale = vscale
     
        val_string = self.infile + '_' + '_'.join(
                [k for k in sorted(vals.values()) if k not in "[]()'"]
                )
        self._info = val_string

        self._set_model(model,merge_vowels)
        self._set_val(
                scale,
                factor,
                restricted_chars,
                pairwise_threshold
                )
        self._make_all_pairs()
        if self.verbose: print("[i] Loaded and calculated all essential values.")
        if score_mode == 'library' and not heuristic:
            try:
                self.scorer = load(open('.'+val_string+'.bin','rb'))
            except:
                self._expand_scorer(
                        runs,
                        modes,
                        ratio
                        )
                dump(self.scorer,open('.'+val_string+'.bin','wb'))
        elif score_mode == 'library':
            try:
                self.scorer = load(open('.'+val_string+'.bin','rb'))
            except:
                self._get_all_pairwise_scores(
                        mode,
                        score_mode='sca',
                        gep_scale = gep_scale,
                        gop = gop,
                        scale = scale,
                        factor = factor
                        )
                # get the cognates for the threshold defined by the heuristic
                self._get_cognates(
                        heuristic
                        )
                if self.verbose: print("[i] Calculated probable cognates.")
                self._expand_scorer(
                        runs,
                        modes,
                        ratio,
                        'sca'
                        )
                dump(self.scorer,open('.'+val_string+'.bin','wb'))

        if self.verbose: print("[i] Created the library.")
        self._get_all_pairwise_scores(
                mode,
                score_mode=score_mode,
                factor = factor,
                gep_scale = gep_scale,
                scale=scale,
                gop = gop
                )
        if self.verbose: print("[i] Calculated pairwise scores.")
        self._get_cognates(threshold,clr_mode)
        if self.verbose: print("[i] Calculated cognates.")
        self._etym_dict()

    def output(
            self,
            fileformat = 'csv',
            filename = None,
            alm_mode = 'lib',
            alm_model = 'sca',
            tree = 'upgma'
            ):
        """
        Write the data to file.

        Parameters
        ----------
        fileformat : { 'lxs', 'star', 'psa', 'msa', 'alm', 'dst' }
            Indicate which data should be written to a file. Select between:

            * 'lxs' -- output in ``lxs``-format,
            * 'star' -- output in ``star``-format (one separate file for each language),
            * 'psa' -- output of all pairwise alignments in ``psa``-format,
            * 'msa' -- output of all cognate sets in separate files in
              ``msa``-format, 
            * 'alm' -- output of all cognate sets in one file in ``alm``-format
              (all data with presumed cognates in aligned format), 
            * 'dst' -- the distance matrix in ``dst``-format (input format for
              distance analyses in the Phylip software package (see
              http://evolution.genetics.washington.edu/phylip/), or
            * 'csv' -- output in ``csv``-format.
        
        filename : str
            Select a specific name for the outfile, otherwise, the name of
            the infile will be taken by default.

        alm_mode : { 'lib', 'prog' }
            Indicate which alignment mode should be used when 'alm' is selected
            as output-format.
        
        """
        # check, if filename is chose as an option
        if not filename:
            filename = self.infile

        outfile = filename + '.' + fileformat

        if fileformat == 'lxs':
            out = open(outfile,'w')
            try:
                self.cogs
                cognates = True
            except:
                cognates = False
            try:
                self.items
                headwords = True
            except:
                headwords = False
            
            if headwords == True:
                if cognates == True:
                    out.write('ID\tWords\t'+'\tCOG\t'.join(self.taxa)+'\tCOG\n')
                    for i in range(len(self.cogs)):
                        out.write(str(self.numbers[i])+'\t'+self.items[i])
                        for j in range(len(self.taxa)):
                            out.write('\t'+self.words[i][j])
                            out.write('\t'+str(self.cogs[i][j]))
                        out.write('\n')
                    out.close()
                elif cognates == False:
                    out.write('ID\tWords\t'+'\t'.join(self.taxa)+'\n')
                    for i in range(len(self.cogs)):
                        out.write(str(self.numbers[i])+'\t'+self.items[i])
                        for j in range(len(self.taxa)):
                            out.write('\t'+self.wordss[i][j])
                        out.write('\n')
                    out.close()
            elif headwords == False:
                if cognates == True:
                    out.write('ID\t'+'\tCOG\t'.join(self.taxa)+'\tCOG\n')
                    for i in range(len(self.cogs)):
                        out.write(str(self.numbers[i]))
                        for j in range(len(self.langs)):
                            out.write('\t'+self.words[i][j])
                            out.write('\t'+str(self.cogs[i][j]))
                        out.write('\n')
                    out.close()
                elif cognates == False:
                    out.write('ID\t'+'\t'.join(self.taxa)+'\n')
                    for i in range(len(self.cogs)):
                        out.write(str(self.numbers[i]))
                        for j in range(len(self.taxa)):
                            out.write('\t'+self.words[i][j])
                        out.write('\n')
                    out.close()
        
        if fileformat == 'star':
            # output the star-files, which can easily be converted into a
            # lexicostatistical list with given cognate judgments (or not)
            try:
                os.mkdir(self.infile+'_star')
                dest = self.infile+'_star/'
            except:
                dest = self.infile+'_star/'
            for i,lang in enumerate(self.taxa):
                #try:
                out = open(dest+lang+'.star','w')
                #except:
                #    out = open(lang+'.star','w')
                for j,word in enumerate(self.words[:,i]):
                    try:
                        cog = str(self.cogs[j,i])
                    except:
                        cog = ''
                    wrd = self.items[j]
                    num = str(self.numbers[j])
                    if word != '-':
                        out.write('\t'.join([num,'---',word,cog])+'\n')
                out.close()

        if fileformat == 'psa':
            # output all pairwise alignments along with the explicit scores of
            # the pairs
            out = open(outfile,'w')
            out.write(self.infile+'\n')
            for i in range(self.width):
                for j in range(self.width):
                    if i < j:

                        for pairA,pairB in self._pairs[i,j]:

                            # mind that the pairs should be the absolute values
                            pairA,pairB = abs(pairA),abs(pairB)

                            # get the alignment
                            almA,almB,scoreAB = self._alignments[pairA,pairB]
                            
                            # get the tokens
                            tokenA = self.tokens[self.idxd[pairA]]
                            tokenB = self.tokens[self.idxd[pairB]]

                            # turn the alignments into ipa strings
                            outA = class2tokens(tokenA,almA)
                            outB = class2tokens(tokenB,almB)

                            # turn all gaps into gap strings in the input data
                            for k,(a,b) in enumerate(zip(almA,almB)):
                                if a == '-':
                                    almA[k] = str(i)+'.X.-'
                                if b == '-':
                                    almB[k] = str(j)+'.X.-'

                            # create a string containing all scores averaged to
                            # round(number,0)
                            scores = [
                                    str(
                                        round(self.scorer[a,b])
                                        ) for a,b in zip(
                                            almA,
                                            almB
                                            )]
                            
                            # get the item
                            item = self[pairA,'i']

                            # get the distance
                            dist = str(self.scores[pairA,pairB])
                            score = str(round(scoreAB))
                            
                            # get the taxa
                            taxA = self._txf.format(self.taxa[i])
                            taxB = self._txf.format(self.taxa[j])

                            # write the data to file
                            out.write(item+'\n')
                            out.write(taxA+'\t'+'\t'.join(outA).encode('utf-8')+'\n')
                            out.write(taxB+'\t'+'\t'.join(outB).encode('utf-8')+'\n')
                            out.write('#'+'\t'+'\t'.join(scores)+'\t'+score+'\n'+dist+'\n\n')
            out.close()

        if fileformat == 'msa':
            try:
                os.mkdir(self.infile+'_msa')
                dest = self.infile+'_msa/'
            except:
                dest = self.infile+'_msa/'

            for key,val in self.etym_dict.items():
                cog = key
                try:
                    entries = ['.'.join(self[i,'t']).encode('utf-8') for i in val[4]]
                except:
                    entries = [self[i] for i in val[4]]
                taxa = [self[i,'l'] for i in val[4]]
                item = self[val[4][0],'i']

                if len(taxa) > 1:
                    mult = _Multiple(entries)
                    mult.lib_align()
                    
                    out = open(dest+self.infile+'_'+str(cog)+'.msa','w')
                    out.write(self.infile+'\n')
                    out.write(item+'\n')
                    for i,alm in enumerate(mult.alm_matrix):
                        out.write(self._txf.format(taxa[i])+'\t')
                        out.write('\t'.join(alm).encode('utf-8')+'\n')
                    out.close()
            
        if fileformat == 'alm':

            out = open(outfile,'w')
            out.write(self.infile+'\n')

            outstring = '{0}\t{1}\t{2}\t{3}\t{4}\n'
            
            # define a previous item for item tracking
            previous_number = 0

            for key,val in sorted(self.etym_dict.items(),key=lambda x:x[0]):
                cog = key
                try:
                    entries = [' '.join(self[i,'t']).encode('utf-8') for i in val[4]]
                except:
                    entries = [self[i] for i in val[4]]
                taxa = [self[i,'l'] for i in val[4]]
                item = self[val[4][0],'i']
                number = self[val[4][0],'N']

                if len(taxa) > 1:
                    mult = _Multiple(entries)
                    if alm_mode == 'lib':
                        mult.lib_align(model=alm_model)
                    else:
                        mult.prog_align(model=alm_model)
                    
                    # redefine previous_number for item tracking
                    if number == previous_number:
                        pass
                    else:
                        previous_number = number
                        out.write('\n')
                    
                    for i,alm in enumerate(mult.alm_matrix):
                        out.write(
                                outstring.format(
                                    key,
                                    self._txf.format(taxa[i]),
                                    item,
                                    number,
                                    '\t'.join(alm).encode('utf-8')
                                    )
                                )
                else:
                    if number == previous_number:
                        pass
                    else:
                        previous_number = number
                        out.write('\n')
                    out.write(
                            outstring.format(
                                key,
                                self._txf.format(taxa[0]),
                                item,
                                number,
                                entries[0],
                                '--'
                                )
                            )
            out.close()
         
        if fileformat == 'alm.patchy':

            out = open(outfile,'w')
            out.write(self.infile+'\n')

            outstring = '{0}.{1}\t{2}\t{3}\t{4}\t{5}\n'
            
            # define a previous item for item tracking
            previous_number = 0

            # create the etymdict
            try:
                self.etym_dict
            except:
                self._etym_dict()

            for key,val in sorted(self.etym_dict.items(),key=lambda x:x[0]):
                
                cog = key
                try:
                    entries = [' '.join(self[i,'t']).encode('utf-8') for i in val[4]]
                except:
                    entries = [self[i] for i in val[4]]
                taxa = [self[i,'l'] for i in val[4]]
                item = self[val[4][0],'i']
                number = self[val[4][0],'N']
                patchies = [self.dataById[i][self._ifs['origin'][0]] for i in
                        val[4]]

                # sort all data according to patchies
                sorter = sorted(
                        zip(range(len(patchies)),patchies),
                        key=lambda x:x[1]
                        )
                sorter = [x[0] for x in sorter]

                entries = [entries[x] for x in sorter]
                taxa = [taxa[x] for x in sorter]
                patchies = [patchies[x] for x in sorter]

                # carry out additional check for double variants with the same
                # cogId (exclude all patchies which are not assigned)
                entries = [entries[x] for x in range(len(entries)) if  patchies[x] != '-']
                taxa = [taxa[x] for x in range(len(taxa)) if patchies[x] != '-']
                patchies = [p for p in patchies if p != '-']

                if len(taxa) > 1:
                    if self.verbose: print(
                            "[i] writing aligned output for COG {0}".format(key)
                            )
                    mult = _Multiple(entries)
                    if alm_mode == 'lib':
                        mult.lib_align(model=alm_model)
                    else:
                        mult.prog_align(model=alm_model)
                    
                    # redefine previous_number for item tracking
                    if number == previous_number:
                        pass
                    else:
                        previous_number = number
                        out.write('\n')
                    
                    for i,alm in enumerate(mult.alm_matrix):
                        out.write(
                                outstring.format(
                                    key,
                                    patchies[i],
                                    self._txf.format(taxa[i]),
                                    item,
                                    number,
                                     '\t'.join(alm).encode('utf-8'),
                                    )
                                )
                #else:
                #    if number == previous_number:
                #        pass
                #    else:
                #        previous_number = number
                #        out.write('\n')
                #    out.write(
                #            outstring.format(
                #                key,
                #                self._txf.format(taxa[0]),
                #                item,
                #                number,
                #                entries[0],
                #                '--',
                #                i
                #                )
                #            )
            out.close()
        
        if fileformat == 'alm2':

            out = open(outfile[:-1],'w')
            out.write(self.infile+'\n')

            outstring = '{0}\t{1}\t{2}\t{3}\t{4}\n'
            
            # define a previous item for item tracking
            previous_number = 0

            for key,val in sorted(self.etym_dict.items(),key=lambda x:x[0]):
                cog = key
                try:
                    entries = ['\t'.join(self[i,'a']).encode('utf-8') for i in val[4]]
                except:
                    entries = [self[i] for i in val[4]]
                taxa = [self[i,'l'] for i in val[4]]
                item = self[val[4][0],'i']
                number = self[val[4][0],'N']

                if len(taxa) > 1:

                    # redefine previous_number for item tracking
                    if number == previous_number:
                        pass
                    else:
                        previous_number = number
                        out.write('\n')
                    
                    for i,alm in enumerate(entries):
                        out.write(
                                outstring.format(
                                    key,
                                    self._txf.format(taxa[i]),
                                    item,
                                    number,
                                    alm 
                                    )
                                )
                else:
                    if number == previous_number:
                        pass
                    else:
                        previous_number = number
                        out.write('\n')
                    out.write(
                            outstring.format(
                                key,
                                self._txf.format(taxa[0]),
                                item,
                                number,
                                entries[0],
                                '--'
                                )
                            )
            out.close()
       
        if fileformat == 'dst':

            out = open(outfile,'w')

            try:
                dist_matrix = self.pairwise_distances()
            except:
                if self.verbose: print("[!] You should conduct cognate judgments first!")
                return

            out.write(' {0}\n'.format(self.width))
            
            for i,taxA in enumerate(self.taxa):
                out.write('{0:<10}'.format(taxA))
                for j,taxB in enumerate(self.taxa):
                    out.write(' {0:.2f}'.format(dist_matrix[i][j]))
                out.write('\n')

            out.close()
        
        # return a neighbor-joining tree
        if fileformat == 'tre':

            out = open(outfile,'w')
            try:
                dist_matrix = self.pairwise_distances()
            except:
                if self.verbose: 
                    print("[!] You should conduct cognate judgments first!")
                return
            if tree == 'upgma':
                tree = upgma(dist_matrix,self.taxa)
            else:
                tree = neighbor(dist_matrix,self.taxa)
            #print(tree)

            out.write(tree)

            out.close()
        
        # return presence-absence-patterns, currently, only straight format
        # without multiple words per cognate set is supported
        if fileformat == 'paps':

            out = open(outfile,'w')
            try:
                self.etym_dict
            except:
                self._etym_dict(include_loans=True)

            # write first line 
            out.write('CogID\t'+'\t'.join(self.taxa)+'\n')
            
            for key in self.etym_dict.keys():

                ids = [x for x in self.etym_dict[key][4]]

                # get languages,ipa, and glossid
                taxa = [self.dataById[x][self._ifs['taxa'][0]] for x in ids]
                words = [self.dataById[x][self._ifs['words'][0]] for x in ids]
                numbers = [self.dataById[x][self._ifs['numbers'][0]] for x in ids]

                # create the cog_id
                cog_id = '{0}.{1}'.format(key,numbers[0])

                # write cog_id to file
                out.write(cog_id)

                # iterate over taxa and write output to file
                for taxon in self.taxa:
                    if taxon in taxa:
                        out.write('\t1')
                    else:
                        out.write('\t0')
                out.write('\n')

            out.close()

        # return etym-dict, currently, only straight format
        # without multiple words per cognate set is supported
        if fileformat == 'pap.ids':

            out = open(outfile,'w')
            try:
                self.etym_dict
            except:
                self._etym_dict(include_loans=True)

            # write first line 
            out.write('CogID\t'+'\t'.join(self.taxa)+'\n')
            
            for key in self.etym_dict.keys():

                ids = [x for x in self.etym_dict[key][4]]

                # get languages,ipa, and glossid
                taxa = [self.dataById[x][self._ifs['taxa'][0]] for x in ids]
                words = [self.dataById[x][self._ifs['words'][0]] for x in ids]
                numbers = [self.dataById[x][self._ifs['numbers'][0]] for x in ids]

                # create the cog_id
                cog_id = '{0}.{1}'.format(key,numbers[0])

                # write cog_id to file
                out.write(cog_id)

                # iterate over taxa and write output to file
                for taxon in self.taxa:
                    if taxon in taxa:
                        idx = ids[taxa.index(taxon)]
                        out.write('\t'+str(idx))
                    else:
                        out.write('\t0')
                out.write('\n')

            out.close()



        if fileformat == 'csv':
            out = open(outfile,'w')
            out.write('# '+self.infile+'\n')
            
            if sum(self.cogs.flatten()) != 0 and 'cogid' not in ','.join(self.dataById[0]).lower():
                out.write('\t'.join(self.dataById[0])+"\tCogID"'\n')
                previous = 0
                for key in sorted(
                        self.idxd.keys(),
                        ):
                    value = self.dataById[key]
                    if value[self._ifs['numbers'][0]] == previous:
                        out.write(
                                '\t'.join(value)+'\t'+str(self.cogs[self.idxd[key]])+'\n'
                                )
                    else:
                        out.write('# \n')
                        out.write(
                                '\t'.join(value)+'\t'+str(self.cogs[self.idxd[key]])+'\n'
                                )
                        previous = value[self._ifs['numbers'][0]]

            else:
                out.write('\t'.join(self.dataById[0])+'\n')
                previous = 0
                for key in sorted(
                        self.idxd.keys(),
                        ):
                    value = self.dataById[key]
                    if value[self._ifs['numbers'][0]] == previous:
                        out.write(
                                '\t'.join(value)+'\n'
                                )
                    else:
                        out.write('# \n')
                        out.write(
                                '\t'.join(value)+'\n'
                                )
                        previous = value[self._ifs['numbers'][0]]
            try:
                out.write('# '+self._info)
            except:
                pass
            out.close()
