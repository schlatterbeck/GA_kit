#!/usr/bin/python3

from math import log
from pga import PGA, PGA_STOP_MAXITER, PGA_STOP_NOCHANGE, PGA_REPORT_STRING, \
                PGA_POPREPL_RTR
from sga import PMBGA

log2 = log (2)

class ECGA (PMBGA) :
    """ Extended Compacat Genetic Algorithm
        The ECGA, originally proposed in [1] replaces the mutation and
        crossover of a genetic algorithm by building and later sampling
        a statistical model of the population in each generation. The
        model uses a minimum description length (MDL) model. A better
        description with an example is found in [2].
        Note that we're using the restricted tournament replacement,
        originally called restricted tournament selection by Harik
        [4,5] and later named restricted tournament replacement by
        Pelikan [6] for population replacement which has recently been
        implemented in the genetic algorithm library PGApack.
        Note that currently we support only bitstring (bool) genes.
        Since we're using model building not the standard GA operations,
        we disable mutation by setting the mutation probability to 0.

        [1] George R. Harik. Linkage learning via probabilistic modeling in
            the ECGA. IlliGAL Report 99010, Illinois Genetic Algorithm
            Lab, January 1999
        [2] Georges R. Harik, Fernando G. Lobo, and Kumara Sastry.
            Linkage learning via probabilistic modeling in the extended
            compact genetic algorithm (ECGA). In Pelikan et al. [3],
            pages 39-61.
        [3] Martin Pelikan, Kumara Sastry, and Erick Cantu-Paz, editors.
            Scalable Optimization via Probabilistic Modelling: From
            Algorithms to Applications, volume 33 of Studies in
            Computational Intelligence. Springer, 2006.
        [4] Georges R. Harik. Finding multiple solutions in problems of
            bounded difficulty. IlliGAL Report 94002, Illinois Genetic
            Algorithm Lab, May 1994.
        [5] Georges R. Harik. Finding multimodal solutions using
            restricted tournament selection. In Eshelman [6], pages
            24-31.
        [6] Larry J. Eshelman, editor. Proceedings of the International
            Conference on Genetic Algorithms (ICGA). Morgan Kaufmann
            July 1995.
        [7] Martin Pelikan. Hierarchical Bayesian Optimization
            Algorithm: Toward a New Generation of Evolutionary
            Algorithms, volume 170 of Studies in Fuzziness and Soft
            Computing. Springer, 2005.
    """

    def build_model (self, p_pop) :
        self.genes = []
        for p in self.parents :
            g = []
            self.genes.append (g)
            for idx in range (len (self)) :
                g.append (self.get_allele (p, p_pop, idx))
        self.partitions = dict (((i,), 1) for i in range (len (self)))
        self.candidates = {}
        for part1 in self.partitions :
            for part2 in self.partitions :
                if part1 == part2 :
                    continue
                self.candidates [part1 + part2] = [part1, part2]
        while self.candidates :
            m = 0
            k = None
            keys = tuple (self.candidates.keys ())
            for c in keys :
                c1, c2 = self.candidates [c]
                d = self.mpm (c1) + self.mpm (c2) - self.mpm (c)
                if d > 0 :
                    if d > m :
                        m = d
                        k = c
                else :
                    self.delete (c)
            if k :
                c1, c2 = self.candidates [k]
                self.delete (c1)
                self.delete (c2)
                keys = tuple (self.candidates.keys ())
                for c in keys :
                    if c1 in self.candidates [c] or c2 in self.candidates [c] :
                        self.delete (c)
                for p in self.partitions :
                    self.candidates [p + k] = [p, k]
                self.partitions [k] = 1
            else :
                assert not self.candidates
    # end def build_model

    def delete (self, partition) :
        if partition in self._mpm_cache :
            del self._mpm_cache [partition]
        if partition in self._probab_cache :
            del self._probab_cache [partition]
        if partition in self.candidates :
            del self.candidates [partition]
        if partition in self.partitions :
            del self.partitions [partition]
    # end def delete

    def entropy (self, partition) :
        d = {}
        for g in self.genes :
            key = []
            for k in partition :
                key.append (g [k])
            key = tuple (key)
            if key not in d :
                d [key] = 0.0
            d [key] += 1
        l = len (self.genes)
        p = []
        for k in d :
            d [k] = d [k] / l
            p.append (d [k])
        self._probab_cache [partition] = d
        return sum (-pk * log (pk) / log2 for pk in p) * len (self.genes)
    # end def entropy

    def mpm (self, partition) :
        if partition not in self._mpm_cache :
            e = self.entropy   (partition)
            s = self.repr_size (partition)
            self._mpm_cache [partition] = e + s
        return self._mpm_cache [partition]
    # end def mpm

    def repr_size (self, partition) :
        return (2.0 ** len (partition) - 1.0) * self.log2n1
    # end def repr_size

    def sample_model (self, c1, c2, c_pop) :
        for i in range (self.pop_size) :
            p = i
            if i == self.pop_size - 1 :
                p = c1
            elif i == self.pop_size - 2 :
                p = c2
            for part in self.partitions :
                r = self.random01 ()
                psum = 0.0
                for k in sorted (self._probab_cache [part].keys ()) :
                    psum += self._probab_cache [part][k]
                    if psum >= r :
                        for idx, bit in zip (part, k) :
                            self.set_allele (p, c_pop, idx, bit) 
                        break
                else :
                    assert False
            #self.print_string (self.file, p, c_pop)
    # end def sample_model

    def post_init (self) :
        self.__super.post_init ()
        self.log2n1  = log (len (self) + 1) / log2
        self._mpm_cache = {}
        self._probab_cache = {}
        self._scache = {}
    # end def post_init

    def print_model (self) :
        for part in self.partitions :
            print (list (part), end = ' ', file = self.file)
        print (file = self.file)
    # end def print_model

# end class ECGA
