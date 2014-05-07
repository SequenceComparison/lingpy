# *-* coding: utf-8 *-*
from data import *
from algorithm.misc import *

class Sequence(object):
    """
    Basic class for handling sound-class sequences.

    Parameters
    ----------

    seq : str
        The input sequence in IPA format.

    model : :py:class:`~lingpy.data.model.Model`
        A :py:class:`~lingpy.data.model.Model` object. Three models are
        predefined and automatically loaded when loading LingPy:

        * 'dolgo' -- a model which is based on the sound-class model of
          :evobib:`Dolgopolsky1986`,

        * 'sca' -- an extension of the "dolgo" sound-class based on
          :evobib:`List2012a`, and

        * 'asjp' -- an independent sound-class model which is based on the
          sound-class model of :evobib:`Brown2008` and the empirical data of
          :evobib:`Brown2011`.

    merge_vowels : bool (default=True)
        Indicate, whether neighboring vowels should be merged into diphtongs,
        or whether they should be kept separated during the analysis.
    
    Attributes
    ----------
    ipa : str
        The original format of the input sequence.

    tokens : list
        A tokenized version of the input sequence.

    classes : str
        A sound-class representation of the input sequence.

    prostring : str
        A string-representation of the prosodic environment of the segments.

    trigram : list
        A list representing the sequence as a trigram.
    
    Examples
    --------
    Initialize a sound-class sequence.

    >>> from lingpy import *
    >>> sca = Sequence('t͡sɔyɡə')
    
    Print out its tokens.
    
    >>> for token in sca.tokens: print(token)
    ... 
    t͡s
    ɔy
    ɡ
    ə
    
    Print out its class-string:

    >>> print(sca.classes)
    CUKE

    Compare the length of the IPA-string with that of the sound-class string.

    >>> len(sca) == len(sca.ipa)
    False

    Access the third element of the sound-class sequence and the IPA-string.

    >>> print(sca[3],sca.ipa[3])
    ə s
    
    Access the prosodic string of the sequence.

    >>> sca.prostring
    '#vC>'

    Access a trigram representation of the sequence.

    >>> sca.trigram
    ['#CU', 'CUK', 'UKE', 'KE$']

    """

    def __init__(
            self,
            seq,
            model = sca,
            merge_vowels = True,
            ):

        # set the string attribute
        self.ipa = seq
        
        # Set the attributes for tokens
        self.tokens = ipa2tokens(
                self.ipa,
                merge_vowels = merge_vowels, 
                diacritics = ipa_diacritics, 
                vowels = ipa_vowels
                )
        
        converter = model.converter

        self.classes = self.tokens[:]
        
        for i in range(len(self.tokens)):
            try:
                self.classes[i] = converter[self.tokens[i]]
            except:
                try:
                    self.classes[i] = converter[self.tokens[i][0]]
                except KeyError:
                    raise ValueError(seq+' cannot be converted'+\
                            'due to <'+self.classes[i]+'>')
        self.classes = ''.join(self.classes)
        
        # retrieve syllabic data from the string
        tokens = self.tokens[:]
        
        for i in range(len(self.tokens)):
            try:
                tokens[i] = art.converter[self.tokens[i]]
            except:
                tokens[i] = art.converter[self.tokens[i][0]]
        
        self.sonar = [int(i) for i in tokens]
        self.prostring = prosodic_string(self.sonar)

        self.trigram = [''.join(i) for i in zip(
            '#'+self.classes[:-1],
            self.classes,
            self.classes[1:]+'$'
            )
            ]

    def __str__(self):

        return '\t'.join(self.tokens)

    def __getitem__(self,x):
        
        try:
            data = x[1]
            x = x[0]
        except:
            data = 't'

        if data == 't':
            return self.tokens[x]
        elif data == 'c':
            return self.classes[x]
        elif data == 'p':
            return self.prostring[x]
        elif data == 'T':
            return self.trigram[x]
        elif data == 's':
            return self.sonar[x]

    def __len__(self):

        return len(self.tokens)


