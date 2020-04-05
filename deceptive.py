from __future__ import print_function
from pga import PGA, PGA_STOP_MAXITER, PGA_REPORT_STRING, PGA_POPREPL_RTR
from pga import PGA_OLDPOP
from rsclib.autosuper import autosuper
from argparse import ArgumentParser
from ecga import ECGA

class Deceptive (ECGA, autosuper) :

    fun = ((3, 1),)
    def __init__ \
        ( self
        , fun         = fun
        , shuffle     = False
        , random_seed = 42
        , popsize     = 1000
        , maxiter     = 1000
        ) :
        self.random_seed = random_seed
        self.fun         = fun
        self.shuffle     = shuffle
        stop_on = [PGA_STOP_MAXITER]
        length  = 0
        for k, n in self.fun :
            length += k * n
        ECGA.__init__ \
            ( self
            , length
            , maximize            = True
            , pop_size            = popsize
            , num_replace         = popsize
            , random_seed         = self.random_seed
            , print_options       = [PGA_REPORT_STRING]
            , stopping_rule_types = stop_on
            , pop_replace_type    = PGA_POPREPL_RTR
            , print_frequency     = 10
            , max_GA_iter         = maxiter
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
        #print (self.funidx)
        self.maxeval = 0.0
        for idxes in self.funidx :
            self.maxeval += len (idxes)
    # end def __init__

    def stop_cond (self) :
        best = self.get_best_index (PGA_OLDPOP)
        eval = self.evaluate (best, PGA_OLDPOP)
        if eval >= self.maxeval :
            return True
        return self.check_stopping_conditions ()
    # end def stop_cond

    def evaluate (self, p, pop) :
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
        return eval
    # end def evaluate

# end class Deceptive

def main () :
    cmd = ArgumentParser ()
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
        ( '-p', '--popsize'
        , type    = int
        , help    = "Population size, default=%(default)s"
        , default = 1000
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

    d = Deceptive \
        ( deceptive_function
        , popsize     = args.popsize
        , random_seed = args.random_seed
        , maxiter     = args.maxiter
        , shuffle     = args.shuffle
        )
    d.run ()
# end def main

if __name__ == '__main__' :
    main ()
