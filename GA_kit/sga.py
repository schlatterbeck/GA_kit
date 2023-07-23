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

import sys
from math import log
from pga  import PGA, PGA_STOP_MAXITER, PGA_STOP_NOCHANGE \
          , PGA_REPORT_STRING, PGA_POPREPL_RTR
from rsclib.autosuper import autosuper

invlog2 = 1.0 / log (2)
def log2 (x) :
    return log (x) * invlog2

class SGA (PGA, autosuper) :
    """ Simple Genetic Algorithm
        Only binary allele are used.
    """

    def __init__ \
        ( self
        , length
        , maximize            = False
        , pop_size            = 100
        , num_replace         = 100 # ignored, always use pop_size
        , random_seed         = 42
        , print_options       = [PGA_REPORT_STRING]
        , stopping_rule_types = [PGA_STOP_NOCHANGE, PGA_STOP_MAXITER]
        , pop_replace_type    = PGA_POPREPL_RTR
        , print_frequency     = 10
        , mutation_prob       = 0.0
        , max_GA_iter         = 1000
        , rtr_window_size     = 0
        , tournament_size     = 2
        , ** kw
        ) :
        if not rtr_window_size :
            rtr_window_size = int (min (pop_size * 0.2, length))
        PGA.__init__ \
            ( self, bool
            , length
            , maximize            = maximize
            , pop_size            = pop_size
            , num_replace         = pop_size
            , random_seed         = random_seed
            , print_options       = print_options
            , stopping_rule_types = stopping_rule_types
            , pop_replace_type    = pop_replace_type
            , print_frequency     = print_frequency
            , max_GA_iter         = max_GA_iter
            , mutation_prob       = mutation_prob
            , crossover_prob      = 1.0
            , rtr_window_size     = rtr_window_size
            , tournament_size     = tournament_size
            )
        self.post_init ()
        self.eval_counter = 0
    # end def __init__

    def post_init (self) :
        pass
    # end def post_init

    def print_string (self, file, p, pop) :
        self.__super.print_string (file, p, pop)
        print ("\nEvaluations: ", self.eval_counter, file = file)
        file.flush ()
    # end def print_string

# end class SGA

class PMBGA (SGA) :
    """ Probabilistic model building GA
        This is a stub, it overrides the mutation to fit the model
        building / sampling into the framework of PGApy.
    """

    def clear_cache (self) :
        pass
    # end def clear_cache

    def post_init (self) :
        self.crossover_count = 0
        self.parents  = []
        self.last_gen = self.get_iteration ()
        self.file     = sys.stdout
    # end def post_init

    def build_model (self, p_pop) :
        self.genes = []
        for p in self.parents :
            g = []
            self.genes.append (g)
            for idx in range (len (self)) :
                g.append (self.get_allele (p, p_pop, idx))
        if getattr (self.__super, 'build_model', None) :
            self.__super.build_model (p_pop)
    # end def build_model

    def crossover (self, p1, p2, p_pop, c1, c2, c_pop) :
        """ Perform crossover from p1, p2 into c1, c2.
            Note that we do not actually perform crossover. We build a
            model from all the parents. Then we sample the model to
            produce the children. This function is called
            self.pop_size / 2 times.
            Note: For this to work the crossover probability must be 1.0
        """
        assert self.crossover_count < self.pop_size
        assert self.get_iteration () == self.last_gen
        self.parents.append (p1)
        self.parents.append (p2)
        self.crossover_count += 2
        if self.crossover_count == self.pop_size :
            assert (self.get_iteration () == self.last_gen)
            print (self.get_iteration ())
            sys.stdout.flush ()
            self.build_model  (p_pop)
            self.sample_model (c1, c2, c_pop)
            self.crossover_count = 0
            self.parents  = []
            self.children = {}
            self.last_gen += 1
            self.clear_cache ()
    # end def crossover

    def print_string (self, file, p, pop) :
        f = self.file
        self.file = file
        self.print_model ()
        self.file.flush ()
        self.__super.print_string (file, p, pop)
        self.file = f
    # end def print_string

    def sample_model (self, c1, c2, c_pop) :
        for i in range (self.pop_size) :
            p = i
            if i == self.pop_size - 1 :
                p = c1
            elif i == self.pop_size - 2 :
                p = c2
            self.sample_individual (p, c_pop)
    # end def sample_model

# end def PMBGA
