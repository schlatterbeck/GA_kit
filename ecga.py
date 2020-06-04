#!/usr/bin/python3
# Copyright (C) 2020 Dr. Ralf Schlatterbeck Open Source Consulting.
# Reichergasse 131, A-3411 Weidling.
# Web: http://www.runtux.com Email: office@runtux.com
# ****************************************************************************
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ****************************************************************************

from __future__ import print_function
from math import log
from sga import PMBGA, log2

def pkey (partition) :
    return (-len (partition), partition [0])

class ECGA (PMBGA) :
    """ Extended Compacat Genetic Algorithm
    """

    def build_model (self, p_pop) :
        self.__super.build_model (p_pop)
        self.partitions = dict (((i,), 1) for i in range (len (self)))
        self.candidates = {}
        for part1 in sorted (self.partitions) :
            for part2 in sorted (self.partitions) :
                if part1 == part2 :
                    continue
                self.candidates [part1 + part2] = [part1, part2]
        while self.candidates :
            m = 0
            k = None
            keys = tuple (self.candidates.keys ())
            for c in sorted (keys) :
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
                for p in sorted (self.partitions) :
                    self.candidates [p + k] = [p, k]
                self.partitions [k] = 1
            else :
                assert not self.candidates
    # end def build_model

    def clear_cache (self) :
        self._mpm_cache    = {}
        self._probab_cache = {}
    # end def clear_cache

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
            for k in sorted (partition) :
                key.append (g [k])
            key = tuple (key)
            if key not in d :
                d [key] = 0.0
            d [key] += 1
        l = len (self.genes)
        p = []
        for k in sorted (d) :
            d [k] = d [k] / l
            p.append (d [k])
        self._probab_cache [partition] = d
        return sum (-pk * log2 (pk) for pk in p) * len (self.genes)
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

    def sample_individual (self, p, pop) :
        for part in sorted (self.partitions) :
            r = self.random01 ()
            psum = 0.0
            for k in sorted (self._probab_cache [part].keys ()) :
                psum += self._probab_cache [part][k]
                if psum >= r :
                    for idx, bit in zip (part, k) :
                        self.set_allele (p, pop, idx, bit)
                    break
            else :
                assert False
        #self.print_string (self.file, p, c_pop)
    # end def sample_model

    def post_init (self) :
        self.__super.post_init ()
        self.log2n1  = log2 (len (self) + 1)
        self.clear_cache ()
    # end def post_init

    def print_model (self) :
        for part in sorted (self.partitions, key = pkey) :
            print (list (part), end = ' ', file = self.file)
        print (file = self.file)
    # end def print_model

# end class ECGA
