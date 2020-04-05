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
            for idx in range (self.string_length) :
                g.append (self.get_allele (p, p_pop, idx))
        self.partitions = dict (((i,), 1) for i in range (self.string_length))
        self.candidates = {}
        for part1 in self.partitions :
            for part2 in self.partitions :
                if part1 == part2 :
                    continue
                self.candidates [part1 + part2] = [part1, part2]
        while self.candidates :
            m = 0
            k = None
            for c in self.candidates.keys () :
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
                for c in self.candidates.keys () :
                    if c1 in self.candidates [c] or c2 in self.candidates [c] :
                        self.delete (c)
                for p in self.partitions :
                    self.candidates [p + k] = [p, k]
                self.partitions [k] = 1
            else :
                assert not self.candidates
        import pdb; pdb.set_trace ()
    # end def build_model

    def delete (self, partition) :
        if partition in self._ecache :
            del self._ecache [partition]
            del self._scache [partition]
        if partition in self.candidates :
            del self.candidates [partition]
        if partition in self.partitions :
            del self.partitions [partition]
    # end def delete

    def mpm (self, partition) :
        e = self.entropy   (partition)
        s = self.repr_size (partition)
        return e + s
    # end def mpm

    def entropy (self, partition) :
        if partition not in self._ecache :
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
            p = [d [k] / l for k in d]
            self._ecache [partition] = \
                sum (-pk * log (pk) / log2 for pk in p) * len (self.genes)
        return self._ecache [partition]
    # end def entropy

    def repr_size (self, partition) :
        if partition not in self._scache :
            self._scache [partition] = \
                (2.0 ** len (partition) - 1.0) * self.log2n1
        return self._scache [partition]
    # end def repr_size

    def post_init (self) :
        self.__super.post_init ()
        self.log2n1  = log (len (self) + 1) / log2
        self._ecache = {}
        self._scache = {}
    # end def post_init
# end class ECGA
