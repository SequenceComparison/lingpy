#! /usr/bin/env python
"""
This is the basic module for the evaluation of automatic analyses.
The module consists of three classes which deal with the evaluation of
automatic analyses (alignments, cognate judgments). The evaluation is based on
the comparison of a gold standard (reference set) with a test set. The
different evaluation measures which can be calculated with help of the
different classes are essentially all based on the calculation of the
proportion to which the test set is similar to the reference set.

The evaluation measures implemented in this module can be divided into two
parts: Those measures which deal with the comparison of automatic alignments,
and those which deal with the comparison of automatic cognate judgments.
"""
from __future__ import division, print_function
from numpy import array,zeros
from ..data import *

class EvalMSA(object):
    """
    Base class for the evaluation of automatic multiple sequence analyses.

    Parameters
    ----------

    gold, test : :py:class:`~lingpy.compare.Multiple`
        The :py:class:`~lingpy.compare.Multiple` objects which shall be
        compared. The first object should be the gold standard and the second
        object should be the test set. 

    Notes
    -----

    Moste of the scores which can be calculated with help of this class are standard
    evaluation scores in evolutionary biology. For a close description on how
    these scores are calculated, see, for example, :evobib:`Thompson1999`,
    :evobib:`List2012`, and :evobib:`Rosenberg2009b.

    See also
    --------
    lingpy.test.evaluate.EvalPSA
    """
    def __init__(self, 
            gold, 
            test
            ):
        self.gold = gold
        self.test = test
    
    def _get_c_scores(self):
        """
        Calculate the c-scores.
        """
        
        almsGold = array(self.gold.alm_matrix).transpose().tolist()
        almsTest = array(self.test.alm_matrix).transpose().tolist()

        commons = len([i for i in almsGold if i in almsTest])
        
        self.cp = commons / len(almsTest)
        self.cr = commons / len(almsGold)
        self.c_ = 2 * commons / (len(almsTest) + len(almsGold))
        try:
            self.cf = 2 * self.cp * self.cr / (self.cp + self.cr)
        except ZeroDivisionError:
            self.cf = 0.0

    def c_score(self,mode=1):
        r"""
        Calculate the column (C) score. 

        Parameters
        ----------

        mode : { 1, 2, 3, 4 }
            Indicate, which mode to compute. Select between:

            1. divide the number of common columns in reference and test
               alignment by the total number of columns in the test alignment
               (the traditional C score described in :evobib:`Thompson1999`,
               also known as "precision" score in applications of information
               retrieval),
            
            2. divide the number of common columns in reference and test
               alignment by the total number of columns in the reference
               alignment (also known as "recall" score in applications of
               information retrieval),
            
            3. divide the number of common columns in reference and test
               alignment by the average number of columns in reference and test
               alignment, or

            4. combine the scores of mode ``1`` and mode ``2`` by computing
               their F-score, using the formula :math:`2 * \frac{pr}{p+r}`,
               where *p* is the precision (mode ``1``) and *r* is the recall
               (mode ``2``).

        Returns
        -------
        score : float
            The C score for reference and test alignments.
        
        Notes
        -----
        The different c-

        See also
        --------
        lingpy.test.evaluate.EvalPSA.c_score

        """
        try:
            self.cp
        except:
            self._get_c_scores()

        if mode == 1:
            return self.cp
        elif mode == 2:
            return self.cr
        elif mode == 3:
            return self.c_
        elif mode == 4:
            return self.cf
        else:
            raise ValueError('The mode you chose is not available.')

    def pir_score(
            self
            ):
        """
        Compute the percentage of identical rows (PIR) score.

        Returns
        -------

        score : float
            The PIR score.

        Notes
        -----
        The PIR score is the number of identical rows (sequences) in reference and test
        alignment divided by the total number of rows.

        See also
        --------
        lingpy.test.evaluate.EvalPSA.pir_score
        """
        almsGold = self.gold.alm_matrix
        almsTest = self.test.alm_matrix
        
        goods = 0.0
        count = 0.0
        for i in range(len(almsGold)):
            if ''.join(almsGold[i]) == ''.join(almsTest[i]):
                goods += 1.0
            
            count += 1.0

        return goods / count

    def sp_score(self,mode=1):
        """
        Calculate the sum-of-pairs (SP) score.

        Parameters
        ----------

        mode : { 1, 2, 3 }
            Indicate, which mode to compute. Select between:

            1. divide the number of common residue pairs in reference and test
               alignment by the total number of residue pairs in the test
               alignment (the traditional SP score described in
               :evobib:`Thompson1999`, also known as "precision" score in
               applications of information retrieval),
            
            2. divide the number of common residue pairs in reference and test
               alignment by the total number of residue pairs in the reference
               alignment (also known as "recall" score in applications of
               information retrieval),
            
            3. divide the number of common residue pairs in reference and test
               alignment by the average number of residue pairs in reference
               and test alignment.

        Returns
        -------

        score : float
            The SP score for gold standard and test alignments.

        Notes
        -----

        The SP score (see :evobib:`Thompson1999`) is calculated by dividing the number of
        identical residue pairs in reference and test alignment by the total
        number of residue pairs in the reference alignment. 

        See also
        --------
        lingpy.test.evaluate.EvalPSA.sp_score
        """
        try:
            return self.sp
        except:
            self._pair_scores()
            return self.sp

    def jc_score(self):
        """
        Calculate the Jaccard (JC) score.

        Returns
        -------
        score : float
            The JC score.

        Notes
        -----
        The Jaccard score (see :evobib:`List2012`) is calculated by dividing the size of
        the intersection of residue pairs in reference and test alignment by
        the size of the union of residue pairs in reference and test alignment.

        See also
        --------
        lingpy.test.evaluate.EvalPSA.jc_score

        """
        
        try:
            return self.jc
        except:
            self._pair_scores()
            return self.jc

    def _pair_scores(
            self,
            weights = False
            ):
        """
        Calculate msa alignment scores by calculating the pairwise scores.
        """
        
        if self.gold == self.test:
            self.sp = 1.0 
            self.o1 = 1.0 
            self.o2 = 1.0 
            self.o_ = 1.0 
            self.jc = 1.0 
            self.cg1 = 1.0 
            self.cg2 = 1.0 
            self.cg_ = 1.0 
            self.cgf = 1.0
            self.pip = 1.0
            return 

        # replace all characters by numbers
        almsGold = zeros(
                (
                    len(
                        self.gold.alm_matrix
                        ),
                    len(
                        self.gold.alm_matrix[0]
                        )
                    )
                )
        almsTest = zeros(
                (
                    len(
                        self.test.alm_matrix
                        ),
                    len(
                        self.test.alm_matrix[0]
                        )
                    )
                )
        
        # select between calculation which is based on an explicit weighting or
        # a calculation which is based on implicit weighting, explicit
        # weighting is done by choosing a specific sound class model and
        # cluster all sequences which are identical, implicit weighting is
        # done otherwise, i.e. identical (pid = 100) sequences are clustered
        # into one sequence in order to avoid getting good scores when there
        # are too many highly identical sequences. 
        if weights:
            self.gold._set_model(weights)
            self._uniseqs = self.gold.int2ext
        
        else:
            self._uniseqs = {}
            for i,seq in enumerate(self.gold.alm_matrix):
                seq = ''.join(seq).replace('-','')
                try:
                    self._uniseqs[seq] += [i]
                except:
                    self._uniseqs[seq] = [i]

        self.weights = {}
        for key in self._uniseqs.keys():
            vals = self._uniseqs[key]
            l = len(vals)
            for val in vals:
                self.weights[val] = (key,1.0 / l)

        # change residues by assining each residue a unique status in both MSAs
        for key in self._uniseqs:
            k = 1
            vals = self._uniseqs[key]
            tmp = []
            for res in self.gold.alm_matrix[vals[0]]:
                if res == '-':
                    tmp.append(0)
                else:
                    tmp.append(k)
                    k += 1
            almsGold[vals[0]] += array(tmp)

        for key in self._uniseqs:
            k = 1
            vals = self._uniseqs[key]
            tmp = []
            for res in self.test.alm_matrix[vals[0]]:
                if res == '-':
                    tmp.append(0)
                else:
                    tmp.append(k)
                    k += 1
            almsTest[vals[0]] += array(tmp)

        # start computation by assigning the variables
        crp = 0.0 # common residue pairs
        trp = 0.0 # residue pairs in test alignment
        rrp = 0.0 # residue pairs in reference alignment
        urp = 0.0 # unique residue pairs in test and reference
        gcrp = 0.0 # common residue pairs including gaps
        gtrp = 0.0 # length of test alignment
        grrp = 0.0 # length of reference alignment
        pip = 0.0 # percentage of identical pairs score

        testL = len(almsTest[0])
        goldL = len(almsGold[0])
        
        # start iteration
        for i,almA in enumerate(almsGold):
            for j,almB in enumerate(almsGold):
                if i < j:
                    gold = zip(almA,almB)
                    test = zip(almsTest[i],almsTest[j])
                    
                    if self.weights[i][0] != self.weights[j][0]:
                        w = self.weights[i][1] * self.weights[j][1]
                    else:
                        w = 0.0

                    # speed up the stuff when sequences are identical
                    if gold == test:
                        tmp = len([x for x in gold if 0 not in x]) * w
                        crp += tmp 
                        trp += tmp 
                        rrp += tmp 
                        urp += tmp 
                        gcrp += testL * w
                        gtrp += testL * w
                        grrp += goldL * w
                        pip += 1 * w
                    
                    else:
                        if [x for x in gold if x != (0,0)] == [y for y in test \
                                if y != (0,0)]:
                            pip += 1 * w

                        crp += len([x for x in gold if x in test and 0 not in
                            x]) * w
                        trp += len([x for x in test if 0 not in x]) * w
                        rrp += len([x for x in gold if 0 not in x]) * w
                        urp += len(set([x for x in gold+test if 0 not in x])) * w
                        gcrp += len([x for x in gold if x in test]) * w
                        gtrp += testL * w
                        grrp += goldL * w
        
        # calculate the scores
        self.sp = crp / rrp
        self.o1 = self.sp
        self.o2 = crp / trp
        self.o_ = 2 * crp / (rrp + trp)
        self.jc = crp / urp
        self.cg1 = gcrp / grrp # recall
        self.cg2 = gcrp / gtrp # precision
        self.cg_ = 2 * gcrp / (grrp + gtrp)
        self.cgf = 2 * (self.cg1 * self.cg2) / (self.cg1 + self.cg2)
        
        l = len(self._uniseqs)
        self.pip = pip / ((l ** 2 - l) / 2)

    def check_swaps(
            self
            ):
        """
        Check for possibly identical swapped sites.

        Returns
        -------

        swap : { -2, -1, 0, 1, 2 }
            Information regarding the identity of swap decisions is coded by
            integers, whereas

            1 -- indicates that swaps are detected in both gold standard and
              testset, whereas a negative value indicates that the positions
              are not identical,

            2 -- indicates that swap decisions are not identical in gold
              standard and testset, whereas a negative value indicates that
              there is a false positive in the testset, and

            0 -- indicates that there are no swaps in the gold standard and the
              testset.
        """
        try:
            swA = self.gold.swap_index
        except:
            swA = False
        try:
            swB = self.test.swap_index
        except:
            swB = False

        if swA and not swB:
            return 2
        elif not swA and swB:
            return -2
        elif swA and swB:
            if swA == swB:
                return 1
            else:# swA != swB:
                return -1
        elif not swA and not swB:
            return 0

class EvalPSA(object):
    """
    Base class for the evaluation of automatic pairwise sequence analyses.

    Parameters
    ----------
    
    gold, test : :py:class:`lingpy.compare.Pairwise`
        The :py:class:`Pairwise <lingpy.compare.Pairwise>` objects which shall be
        compared. The first object should be the gold standard and the second
        object should be the test set.    
    
    Notes
    -----

    Moste of the scores which can be calculated with help of this class are standard
    evaluation scores in evolutionary biology. For a close description on how
    these scores are calculated, see, for example, :evobib:`Thompson1999`,
    :evobib:`List2012`, and :evobib:`Rosenberg2009b`.

    See also
    --------
    lingpy.test.evaluate.EvalMSA
    """
    def __init__(self,gold,test):

        self.gold = gold
        self.test = test

    def pir_score(
            self,
            mode = 1
            ):
        """
        Compute the percentage of identical rows (PIR) score.

        Parameters
        ----------

        mode : { 1, 2 }
            Select between mode ``1``, where all sequences are compared with
            each other, and mode ``2``, where only whole alignments are
            compared.  

        Returns
        -------

        score : float
            The PIR score.

        Notes
        -----
        The PIR score is the number of identical rows (sequences) in reference and test
        alignment divided by the total number of rows.

        See also
        --------
        lingpy.test.evaluate.EvalMSA.pir_score
        """
        
        score = 0.0
        count = 0.0

        if mode == 1:
            for i,alms in enumerate(self.gold.alignments):
                if self.test.alignments[i][0] == alms[0]:
                    score += 0.5
                if self.test.alignments[i][1] == alms[1]:
                    score += 0.5
                count += 1.0
        elif mode == 2:
            for i,alms in enumerate(self.gold.alignments):
                tmp = 0
                if self.test.alignments[i][0] == alms[0]:
                    tmp = 1
                if self.test.alignments[i][1] == alms[1]:
                    tmp += 1
                if tmp == 2:
                    score += 1.0
                count += 1.0

        return score / count

    def _pairwise_column_scores(
            self
            ):
        """
        Compute the different column scores for pairwise alignments. The method
        returns the precision, the recall score, and the f-score, following the
        proposal of Bergsma and Kondrak (2007), and the column score proposed
        by Thompson et al. (1999).
        """

        # the variables which store the different counts
        
        crp = 0.0 # number of common residue pairs in reference and test alm.
        rrp = 0.0 # number of residue pairs in reference alignment 
        trp = 0.0 # number of residue pairs in test alignment
        urp = 0.0 # number of unique residue pairs in reference and test alm.
        
        gtrp = 0.0 # number of residue pairs (including gaps) in test alm.
        grrp = 0.0 # number of residue pairs (including gaps) in reference alm.
        gcrp = 0.0 # number of common residue pairs (including gaps) in r and t

        self.sps_list = []
        self.cs_list = []

        for i,alms in enumerate(self.gold.alignments):
            zipsA = zip(
                    self.gold.alignments[i][0],
                    self.gold.alignments[i][1]
                    )
            zipsB = zip(
                    self.test.alignments[i][0],
                    self.test.alignments[i][1]
                    )
            
            # replace all residues in reference and test alignment with ids
            pairsGold = []
            j,k = 1,1
            for a,b in zipsA:
                x,y = 0,0
                if a != '-':
                    x = j
                    j += 1
                if b != '-':
                    y = k
                    k += 1
                pairsGold.append((x,y))

            pairsTest = []
            j,k = 1,1
            for a,b in zipsB:
                x,y = 0,0
                if a != '-':
                    x = j
                    j += 1
                if b != '-':
                    y = k
                    k += 1
                pairsTest.append((x,y))

            # calculate the number of residues in crp, rrp, and trp
            commons = len([x for x in pairsTest if x in pairsGold and 0 not in x])
            nogaps = len([x for x in pairsGold if 0 not in x])
            crp += len([x for x in pairsTest if x in pairsGold and 0 not in x])
            rrp += len([x for x in pairsGold if 0 not in x])
            trp += len([x for x in pairsTest if 0 not in x])
            urp += len(set([x for x in pairsGold+pairsTest if 0 not in x]))

            grrp += len(pairsGold)
            gtrp += len(pairsTest)
            gcrp += len([x for x in pairsTest if x in pairsGold])
            
            # fill in list with exact scores
            commons = len([x for x in pairsTest if x in pairsGold and 0 not in x])
            nogaps = len([x for x in pairsGold if 0 not in x])
            nogaps2 = nogaps + len([x for x in pairsTest if 0 not in x])
            columns = len([x for x in pairsTest if x in pairsGold])

            if nogaps != 0:
                self.sps_list.append((2 * commons) / (nogaps + nogaps2)) #new
            elif nogaps == 0 and commons == 0:
                self.sps_list.append(1)
            else:
                self.sps_list.append(0)
            self.cs_list.append((2 * columns)/(len(pairsGold)+len(pairsTest)))#new
        
        # calculate the scores        
        self.sop = crp / rrp
        self.jac = crp / urp
        self.o1 = self.sop
        self.o2 = crp / trp
        self.o_ = 2 * crp / (rrp + trp)
        self.precision = gcrp / gtrp
        self.recall = gcrp / grrp
        self.pic = self.precision
        self.fscore = 2 * (self.precision * self.recall) / (self.precision + \
                self.recall)

    def c_score(
            self
            ):
        """
        Calculate column (C) score. 

        Returns
        -------
        score : float
            The C score for reference and test alignments.
        
        Notes
        -----
        The C score, as it is described in :evobib:`Thompson1999`, is calculated by
        dividing the number of columns which are identical in the gold
        standarad and the test alignment by the total number of columns in the
        test alignment.

        See also
        --------
        lingpy.test.evaluate.EvalMSA.c_score

        """
        try:
            return self.pic
        except:
            self._pairwise_column_scores()
            return self.pic

    def sp_score(
            self
            ):
        """
        Calculate the sum-of-pairs (SP) score.

        Returns
        -------

        score : float
            The SP score for reference and test alignments.

        Note
        ----
        The SP score (see :evobib:`Thompson1999`) is calculated by dividing the number of
        identical residue pairs in reference and test alignment by the total
        number of residue pairs in the reference alignment. 
        
        See also
        --------
        lingpy.test.evaluate.EvalMSA.sp_score

        """
        try:
            return self.sop
        except:
            self._pairwise_column_scores()
            return self.sop

    def jc_score(
            self
            ):
        """
        Calculate the Jaccard (JC) score.

        Returns
        -------
        score : float
            The JC score.

        Notes
        -----

        The Jaccard score (see :evobib:`List2012`) is calculated by dividing the size of
        the intersection of residue pairs in reference and test alignment by
        the size of the union of residue pairs in reference and test alignment.
        
        See also
        --------
        lingpy.test.evaluate.EvalMSA.jc_score

        """
        try:
            return self.jac
        except:
            self._pairwise_column_scores()
            return self.jac

    def diff(
            self,
            filename = None
            ):
        """
        Write all differences between two sets to a file.

        Parameters
        ----------

        filename : str (default='eval_psa_diff')
            Default

        """
        if filename:
            out = open(filename+'_diff.psa','w')
        else:
            out = open('eval_psa_diff.psa','w')

        for i,(a,b) in enumerate(zip(self.gold.alignments,self.test.alignments)):
            
            g1,g2,g3 = a
            t1,t2,t3 = b
            maxL = 1 + max([len(g1),len(t1)])
            if g1 != t1 or g2 != t2:
                taxA,taxB = self.gold.taxa[i]
                seq_id = self.gold.seq_ids[i]
                out.write('{0}\n{1}\t{2}\n{3}\t{4}\n{5}\n{1}\t{6}\n{3}\t{7}\n\n'.format(
                    seq_id,
                    taxA,
                    '\t'.join(g1).encode('utf-8'),
                    taxB,
                    '\t'.join(g2).encode('utf-8'),
                    '\t'.join(['==' for x in range(maxL)]),
                    '\t'.join(t1).encode('utf-8'),
                    '\t'.join(t2).encode('utf-8'),
                    ))
        out.close()

class EvalLXS(object):
    """
    Basic class for comparing automatic lexicostatistical analyses.

    Parameters
    ----------

    gold, test : :py:class:`lingpy.lexstat.LexStat`
        The :py:class:`LexStat <lingpy.lexstat.LexStat>` objects which shall be
        compared. The first object should be the gold standard and the second
        object should be the test set.    
    
    """

    def __init__(self,gold,test,include_loans=False):

        self.gold = gold
        self.test = test
        self.include_loans = include_loans

    def pid_score(self):
        """
        Compute the pairwise identical decisions (PID) score.

        Returns
        -------
        score : float
            The PID score for reference and test set.

        Notes
        -----
        The PID score (see :evobib:`List2012b`), is calculated by dividing the number of
        identical pairwise decisions in reference and test set by the total
        number of pairwise decisions.

        See also
        --------
        lingpy.test.evaluate.EvalLXS.pind_score
        lingpy.test.evaluate.EvalLXS.pipd_score


        """
        try:
            return self.pid
        except:
            self._compare_cognate_pairs()
            return self.pid

    def pind_score(self):
        """
        Compute the pairwise identical negative decisions (PIND) score.

        Returns
        -------
        score : float
            The PIND score for reference and test set.

        Notes
        -----
        The PIND score (see :evobib:`List2012b`), is calculated by dividing the
        number of identical pairwise negative decisions in reference and test
        set by the total number of pairwise negative decisions in the reference
        set.

        See also
        --------
        lingpy.test.evaluate.EvalLXS.pid_score
        lingpy.test.evaluate.EvalLXS.pipd_score

        """
        try:
            return self.pind
        except:
            self._compare_cognate_pairs()
            return self.pind

    def pipd_score(self):
        """
        Compute the pairwise identical positive decisions (PIPD) score.

        Returns
        -------
        score : float
            The PIPD score for reference and test set.

        Notes
        -----
        The PIPD score (see :evobib:`List2012b`), is calculated by dividing the number of
        identical pairwise positive decisions in reference and test set by the total
        number of pairwise positive decisions in the reference set.

        See also
        --------
        lingpy.test.evaluate.EvalLXS.pid_score
        lingpy.test.evaluate.EvalLXS.pind_score

        """
        try:
            return self.pipd
        except:
            self._compare_cognate_pairs()
            return self.pipd

    def _compare_cognate_pairs(self):
        
        # create the dictionaries of cognates
        # pass pairs of test to gold standard
        self.gold._pairs = self.test._pairs
    
        gold_standard = self.gold._make_cognate_pairs(
                include_loans=self.include_loans
                )
        test = self.test._make_cognate_pairs(
                include_loans=self.include_loans
                )
    
        general_score = 0.0
        counter = 0
        truep = 0.0
        falsp = 0.0
        truew = 0.0
        falsw = 0.0
    
        for key in gold_standard:
            a = gold_standard[key]
            b = test[key]
            if a == b:
                general_score += 1.0
                counter += 1
            else:
                counter += 1
    
            if a == 1 and b == 1:
                truep += 1.0
            elif a == 1 and b == 0:
                falsw += 1.0
            elif a == 0 and b == 1:
                falsp += 1.0
            elif a == 0 and b == 0:
                truew += 1.0
    
        # calculate general score, true and false positives, etc.
        pipd = truep / (truep + falsw)
        pind = truep / (truep + falsp)

        self.ID = (truep + truew) / (truep + truew + falsp + falsw)
        
        pid = general_score / counter
    
        false_positives = falsp / counter
        false_negatives = falsw / counter
        
        self.pid = 2 * (pipd * pind) / (pipd + pind) #pid
        self.pipd = pipd
        self.pind = pind
        self.fp = false_positives
        self.fn = false_negatives
        self.number_of_pairs = counter
        self.number_of_trues = truep+falsw
        self.number_of_wrongs = truew+falsp

        # precision and recall can be easily defined from the number of cognate
        # pairs, etc.
        self._pair_precision = pind
        self._pair_recall = pipd
        self._pair_fscore = self.pid

    def pair_precision(self):
        try:
            return self._pair_precision
        except:
            self._set_scores()
            return self._pair_precision

    def pair_recall(self):
        try:
            return self._pair_recall
        except:
            self._set_scores()
            return self._pair_recall
    def pair_fscore(self):
        try:
            return self._pair_fscore
        except:
            self._set_scores()
            return self._pair_fscore
    
    def set_precision(self):
        """
        Compute the set precision (SPR).

        Returns
        -------
        score : float
            The SPR for reference and test set.

        Notes
        -----
        The set precision (see :evobib:`Bergsma2007`) is defined as the number of identical
        cognate sets in reference and test set divided by the total number of
        cognate sets in the test set.

        See also
        --------
        lingpy.test.evaluate.EvalLXS.set_recall
        lingpy.test.evaluate.EvalLXS.set_fscore

        """
        try:
            return self.pre
        except:
            self._cognate_set_scores()
            return self.pre

    def set_recall(self):
        """
        Compute the set recall (SRE).

        Returns
        -------
        score : float
            The SRE for reference and test set.

        Notes
        -----
        The set recall (see :evobib:`Bergsma2007`) is defined as the number of identical
        cognate sets in reference and test set divided by the total number of
        cognate sets in the reference set.

        See also
        --------
        lingpy.test.evaluate.EvalLXS.set_precision
        lingpy.test.evaluate.EvalLXS.set_fscore

        """
        try:
            return self.rec
        except:
            self._cognate_set_scores()
            return self.rec
    
    def set_fscore(self):
        r"""
        Compute the set f-score (FSC).

        Returns
        -------
        score : float
            The FSC for reference and test set.

        Notes
        -----
        The set f-score (see :evobib:`Bergsma2007`) is calculated with help of the formula

        .. math:: 

            2 \frac{pr}{p+r},

        where `p` is the lGold set precision and `r` is the set recall.

        See also
        --------
        lingpy.test.evaluate.EvalLXS.set_precision
        lingpy.test.evaluate.EvalLXS.set_recall
        """
        try:
            return self.fsc
        except:
            self._cognate_set_scores()
            return self.fsc

    def _set_scores(self):
        """
        Calculate set precision following Bergsma and Kondrak (2007).
        """
        try:
            edGold = self.gold.etym_dict
        except:
            self.gold._etym_dict(
                    include_loans=self.include_loans
                    )
            edGold = self.gold.etym_dict
        edTest = self.test.etym_dict
    
        new_edGold,new_edTest = [],[]
    
        for key in edGold:
            new_edGold.append(
                    set(
                        edGold[key][4]
                        )
                    )
        for key in edTest:
            new_edTest.append(
                    set(
                        edTest[key][4]
                        )
                    )
    
        lGold = len(new_edGold)
        lTest = len(new_edTest)
        
        # calculate the number of sets in the gold standard which are also present
        # in the test set
        intersectionTestGold = len([i for i in new_edGold if i in new_edTest])
        
        # calculate precision and recall 
        precision = intersectionTestGold / lTest
        recall = intersectionTestGold / lGold
    
        fscore =  2 * ((precision * recall) / (precision + recall))

        self.precision = precision
        self.recall = recall
        self.fscore = fscore
        self.common_sets = intersectionTestGold
        self.gold_sets = lGold
        self.test_sets = lTest

    def _cognate_set_scores(self):
        """
        Calculate relaxed form of precision and recall, where not set identity
        but set proportion is calculated.
        """
        
        # get the etym_dicts
        self.test._etym_dict(include_loans=self.include_loans)
        self.gold._etym_dict(include_loans=self.include_loans)

        # get sets in gold and test
        lGold = [self.gold.etym_dict[k][4] for k in sorted(self.gold.etym_dict.keys())]
        lTest = [self.test.etym_dict[k][4] for k in sorted(self.test.etym_dict.keys())]

        # start comparing from gold perspective (equals the recall)
        recall = []
        bcr = [] # bcube recall
        for g in lGold:
            if len(g) > 1:
                # get the cognate ids in the testset
                l = [self.test[i,'C'] for i in g]

                if g in lTest:#len(set(l)) == 1:
                    score = 1.0
                else:
                    score = 0.0
                recall.append(score)

                # calculate bcube recall
                # get the length of all clusters
                L = len(l)
                # get the recall
                for word in l:
                    bcr.append(l.count(word) / L)
            else:
                if g in lTest:
                    score = 1.0
                else:
                    score = 0.0
                recall.append(score)
                bcr.append(1.0)

        # start comparing from gold perspective (equals the recall)
        precision = []
        bcp = [] # bcube precision
        for g in lTest:
            if len(g) > 1:
                
                l = [self.gold[i,'C'] for i in g]

                if g in lGold: #len(set(l)) == 1:
                    score = 1.0
                else:
                    score = 0.0

                precision.append(score)
                
                # calculate the bcube precision following Hauer & Kondrak 2011
                # get the length of the cluster
                L = len(l)
                # iterate over each word in l and append its item-precision
                for word in l:
                    # calculate number of "real" cognates
                    bcp.append(l.count(word) / L)
            else:
                if g in lGold:
                    score = 1.0
                else:
                    score = 0.0
                bcp.append(1.0)
                precision.append(score)

        
        self.rec = len([i for i in lTest if i in lGold]) / len(lGold) #array(recall).mean()
        self.pre = len([i for i in lTest if i in lGold]) / len(lTest) #array(precision).mean()
        self.css = 2 * ((self.pre * self.rec) / (self.pre + self.rec))
        self.fsc = self.css

        self.bcp = array(bcp).mean()
        self.bcr = array(bcr).mean()
        self.bcf = 2 * ((self.bcp * self.bcr) / (self.bcp + self.bcr))
    
    def precision(self):
        try:
            return self.bcp
        except:
            self._cognate_set_scores()
            return self.bcp

    def recall(self):
        try:
            return self.bcr
        except:
            self._cognate_set_scores()
            return self.bcr

    def cs_score(self):
        try:
            return self.bcf
        except:
            self._cognate_set_scores()
            return self.bcf

    def bcubed_precision(self):
        try:
            return self.bcp
        except:
            self._cognate_set_scores()
            return self.bcp

    def bcubed_recall(self):
        try:
            return self.bcr
        except:
            self._cognate_set_scores()
            return self.bcr

    def bcubed_fscore(self):
        try:
            return self.bcf
        except:
            self._cognate_set_scores()
            return self.bcf
  
    def compare_pairwise_decisions(self,threshold):
        """
        Compute the number of identical decisions for pairwise scores.

        Parameters
        ----------

        threshold : float
            The threshold which determines the cognacy decisions in the test
            set.

        Returns
        -------

        scores : tuple
            A tule containing the scores for true positives, true negatives and
            the general percentage of identical decisions (PID) score.
        """
        gold_standard = self.gold._make_cognate_pairs(self.include_loans)
        test = self.test.scores
        
        general_score = 0.0
        counter = 0.0
        truep = 0.0
        falsp = 0.0
        truew = 0.0
        falsw = 0.0
    
        for key in gold_standard:
            a = gold_standard[key]
            b = test[key]
            if a == 1 and b <= threshold:
                general_score += 1.0
                counter += 1
            elif a == 0 and b > threshold:
                general_score += 1.0
                counter += 1
            else:
                counter += 1
    
            if a == 1 and b <= threshold:
                truep += 1.0
            elif a == 1 and b > threshold:
                falsw += 1.0
            elif a == 0 and b <= threshold:
                falsp += 1.0
            elif a == 0 and b > threshold:
                truew += 1.0
    
        # calculate general score, true and false positives, etc.
        true_cogs = truep / (truep + falsw)
        true_wrgs = truew / (truew + falsp)
        general_score = general_score / counter
        false_positives = falsp / counter
        false_negatives = falsw / counter
    
        return true_cogs,true_wrgs,general_score

