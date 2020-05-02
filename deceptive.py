#!/usr/bin/python3
from __future__ import print_function
from pga import PGA, PGA_STOP_MAXITER, PGA_REPORT_STRING, PGA_POPREPL_RTR
from pga import PGA_OLDPOP
from rsclib.autosuper import autosuper
from argparse import ArgumentParser
from sga  import SGA
from ecga import ECGA
from hboa import HBOA

class Deceptive (autosuper) :

    fun = ((3, 1),)
    def __init__ \
        ( self
        , fun             = fun
        , shuffle         = False
        , random_seed     = 42
        , popsize         = 1000
        , maxiter         = 1000
        , tournament_size = 2
        , rtr_window_size = 0
        , s_penalty       = 2.0
        , min_split       = 0
        , max_parent      = 0
        ) :
        self.random_seed     = random_seed
        self.fun             = fun
        self.shuffle         = shuffle
        self.s_penalty       = s_penalty
        self.min_split       = min_split
        self.max_parent      = max_parent
        self.tournament_size = tournament_size
        self.rtr_window_size = rtr_window_size
        stop_on = [PGA_STOP_MAXITER]
        length  = 0
        for k, n in self.fun :
            length += k * n
        self.__super.__init__ \
            ( length
            , maximize            = True
            , pop_size            = popsize
            , num_replace         = popsize
            , random_seed         = self.random_seed
            , print_options       = [PGA_REPORT_STRING]
            , stopping_rule_types = stop_on
            , pop_replace_type    = PGA_POPREPL_RTR
            , print_frequency     = 10
            , max_GA_iter         = maxiter
            , tournament_size     = tournament_size
            , rtr_window_size     = rtr_window_size
            )

        indexes = list (range (len (self)))
        self.funidx = []
        if self.shuffle :
            l = len (self)
            for k in range (l) :
                j = self.random_interval (0, l - 1)
                indexes [k], indexes [j] = indexes [j], indexes [k]
        idx = 0
        for bits, nfunc in self.fun :
            for n in range (nfunc) :
                a = []
                self.funidx.append (a)
                for b in range (bits) :
                    a.append (indexes [idx])
                    idx += 1
        print ("Optimizing:")
        print ("Random seed:     %s" % random_seed)
        print ("Population size: %s" % popsize)
        print ("Tournament size: %s" % tournament_size)
        print ("RTR window size: %s" % rtr_window_size)
        if isinstance (self, HBOA) :
            print ("HBOA S-Penalty:  %s" % s_penalty)
            print ("HBOA Min-Split:  %s" % min_split)
        print ("Functions:", self.funidx)
        self.maxeval = 0.0
        for idxes in self.funidx :
            self.maxeval += len (idxes)
    # end def __init__

    def stop_cond (self) :
        best = self.get_best_index (PGA_OLDPOP)
        eval = self.evaluate (best, PGA_OLDPOP, count_eval = False)
        if eval >= self.maxeval :
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

    def evaluate (self, p, pop, count_eval = True) :
        eval = 0.0
        for indexes in self.funidx :
            l = len (indexes)
            v = 0.0
            for idx in indexes :
                a = self.get_allele (p, pop, idx)
                if a :
                    v += 1
            if v == 0 :
                eval += l
            else :
                eval += v - 1.0
        if count_eval :
            self.eval_counter += 1
        return eval
    # end def evaluate

# end class Deceptive

class Dec_SGA (Deceptive, SGA) :
    pass

class Dec_ECGA (Deceptive, ECGA) :
    pass

class Dec_HBOA (Deceptive, HBOA) :
    pass

def main () :
    cmd = ArgumentParser ()
    classes = ('SGA', 'ECGA', 'HBOA')
    cmd.add_argument \
        ( '-c', '--class'
        , dest    = 'cls'
        , help    = "Class to use, one of %s, default=%%(default)s" % (classes,)
        , default = 'SGA'
        )
    cmd.add_argument \
        ( '-d', '--deceptive-function'
        , help    = "Add deceptive function with length/count"
        , action  = "append"
        )
    cmd.add_argument \
        ( '-m', '--maxiter'
        , type    = int
        , help    = "Maximum number of generations, default=%(default)s"
        , default = 1000
        )
    cmd.add_argument \
        ( '--max-parent'
        , type    = int
        , help    = "Maximum number of parents or 0 for no maximum, "
                    "default:%(default)s, only used for HBOA"
        , default = 0
        )
    cmd.add_argument \
        ( '--min-split'
        , type    = int
        , help    = "Minimum number of samples for split, default:%(default)s, "
                    "only used for HBOA"
        , default = 0
        )
    cmd.add_argument \
        ( '-p', '--popsize'
        , type    = int
        , help    = "Population size, default=%(default)s"
        , default = 1000
        )
    cmd.add_argument \
        ( '--rtr-window-size'
        , type    = int
        , help    = "Window-size for RTR"
        , default = 0
        )
    cmd.add_argument \
        ( '--s-penalty'
        , type    = float
        , help    = "Penalty to apply to the cutoff, default: tournament-size, "
                    "only used for HBOA"
        , default = 0.0
        )
    cmd.add_argument \
        ( '--tournament-size'
        , type    = int
        , help    = "Number of participant in tournament selection, "
                    "default= %(default)s"
        , default = 2
        )
    cmd.add_argument \
        ( '-s', '--shuffle'
        , help    = "Shuffle genes of deceptive functions"
        , action  = "store_true"
        )
    cmd.add_argument \
        ( '-R', '--random-seed'
        , type    = int
        , help    = "Random number seed, default=%(default)s"
        , default = 42
        )
    args = cmd.parse_args ()
    if args.deceptive_function :
        deceptive_function = []
        for d in args.deceptive_function :
            a, b = (int (i) for i in d.split ('/'))
            deceptive_function.append ((a, b))
    else :
        deceptive_function = ((5, 20),)

    cls = globals () ['Dec_' + args.cls]
    if not args.s_penalty :
        args.s_penalty = float (args.tournament_size)

    d = cls \
        ( deceptive_function
        , popsize         = args.popsize
        , random_seed     = args.random_seed
        , maxiter         = args.maxiter
        , shuffle         = args.shuffle
        , rtr_window_size = args.rtr_window_size
        , tournament_size = args.tournament_size
        , s_penalty       = args.s_penalty
        , min_split       = args.min_split
        , max_parent      = args.max_parent
        )
    d.run ()
# end def main

if __name__ == '__main__' :
    main ()
