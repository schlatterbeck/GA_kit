#!/usr/bin/python3

from pga import PGA, PGA_STOP_MAXITER, PGA_STOP_NOCHANGE, PGA_REPORT_STRING, \
                PGA_POPREPL_RTR
from rsclib.autosuper import autosuper

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
        ) :
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
            )
        self.string_length   = length
        self.post_init ()
    # end def __init__

    def post_init (self) :
        pass
    # end def post_init

# end class SGA

class PMBGA (SGA) :
    """ Probabilistic model building GA
        This is a stub, it overrides the mutation to fit the model
        building / sampling into the framework of PGApy.
    """

    def post_init (self) :
        self.crossover_count = 0
        self.parents  = []
        self.children = []
        self.last_gen = self.get_iteration ()
    # end def post_init

    def crossover (self, p1, p2, p_pop, c1, c2, c_pop) :
        """ Perform crossover from p1, p2 into c1, c2.
            Note that we do not actually perform crossover. We build a
            model from all the parents and store the children indexes.
            Then we sample the model to produce the children. This
            function is called self.pop_size / 2 times.
        """
        assert self.crossover_count < len (self)
        self.children.append (c1)
        self.children.append (c2)
        self.parents.append (p1)
        self.parents.append (p2)
        self.crossover_count += 2
        if self.crossover_count == len (self) :
            self.build_model  (p_pop)
            self.sample_model (c_pop)
            self.crossover_count = 0
            self.parents  = []
            self.children = []
            self.last_gen = self.get_iteration ()
    # end def crossover

# end def PMBGA
