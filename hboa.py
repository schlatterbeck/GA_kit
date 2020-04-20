#!/usr/bin/python3

from sys  import stderr
from math import log, lgamma
from sga  import PMBGA, log2

class Bayesian_Network (object) :

    def __init__ \
        ( self, hboa, genes
        , lgamma    = lgamma
        , do_debug  = 0
        , verbose   = 0
        , maxparent = 0
        , min_n     = 30 # do not split below that number
        ) :
        self.hboa      = hboa
        self.nodecount = len (genes [0])
        self.n         = len (genes)
        self.genes     = genes
        self.nodes     = {}
        self.roots     = {}
        self.cutoff    = log2 (self.n) / 2.0
        self.cutoff    = log2 (self.n)
        self.do_debug  = do_debug
        self.verbose   = verbose
        self.lgamma    = lgamma
        self.maxparent = maxparent
        self.min_n     = min_n
        for bitidx in range (self.nodecount) :
            node = BNode (self, bitidx)
            self.nodes [node] = 1
            self.roots [node] = 1
        # The initial candidates already need the full list of nodes
        for n in self.nodes :
            n.add_initial_candidates ()
        nsplit = 0
        while 1 :
            maxgain = -1
            leave   = None
            cidx    = -1
            for node in self.nodes :
                for l in node.leaves () :
                    keys = list (l.candidates.keys ())
                    #for c in l.candidate_iter () :
                    for k in keys :
                        c = l.candidates [k]
                        assert c.idx == k
                        if c.feasible () :
                            if c.gain > maxgain :
                                maxgain = c.gain
                                leave   = l
                                cidx    = c.idx
                            # Since candidates are returned best-first
                            # we can stop after the first found
                            #break
                        else :
                            del l.candidates [k]
                            #l.del_candidate (c)
            if maxgain <= 0 :
                break
            if self.verbose :
                print ("maxgain: %2.2f" % maxgain, end = ' ')
                print \
                    ( "%4.2f -> %4.2f %4.2f"
                    % ( leave.score
                      , leave.candidates [cidx].children [0].score
                      , leave.candidates [cidx].children [1].score
                      )
                    )
            nsplit += 1
            leave.split (cidx)
        print ("nsplit: %s" % nsplit)
    # end def __init__

    def debug (self, *args, **kw) :
        if self.do_debug :
            print (*args, **kw, file = stderr)
    # end def debug

    def _sample_node (self, node, d) :
        dn = node.dnode
        while not isinstance (dn, DLeaf) :
            dn = dn.children [d [dn.idx]]
        assert dn.idx == node.idx
        bit = self.hboa.random_flip (dn.p)
        d [dn.idx] = bit
        for c in node.children :
            if c.idx in d :
                continue
            for parent in c.parents :
                if parent.idx not in d :
                    break
            else :
                self._sample_node (c, d)
    # end def _sample_node

    def sample_model (self) :
        d = {}
        for r in self.roots :
            self._sample_node (r, d)
        return d
    # end def sample_model

# end class Bayesian_Network

class BNode (object) :

    def __init__ (self, net, idx) :
        self.net      = net
        self.debug    = net.debug
        self.verbose  = net.verbose
        self.genes    = net.genes
        self.idx      = idx
        self.parents  = {}
        self.children = {}
        self.lvl      = 0
        self.dnode    = DLeaf (self, self, 0, self.genes)
        self.rank     = 0
    # end def __init__

    def add_initial_candidates (self) :
        for node in self.net.nodes :
            if node.idx == self.idx :
                continue
            self.dnode.try_add_candidate (node)
    # end def add_initial_candidates

    def append_parent (self, node) :
        self.debug ("append_parent: %d %d" % (self.idx, node.idx))
        assert node not in self.children
        assert self not in node.parents
        self.parents  [node] = 1
        node.children [self] = 1
        if self in self.net.roots :
            del self.net.roots [self]
    # end def append_parent

    def leaves (self, root = None) :
        """ The root can also be a leaf
        """
        if not root :
            root = self.dnode
        if isinstance (root, DLeaf) :
            yield (root)
        else :
            for c in root.children :
                if isinstance (c, DLeaf) :
                    yield (c)
                else :
                    for cx in self.leaves (c) :
                        yield (cx)
    # end def leaves

    def is_transitive_parent (self, node) :
        for c in self.children :
            if node is c :
                return True
            if c.is_transitive_parent (node) :
                return True
        return False
    # end def is_transitive_parent

    def may_append_parent (self, node) :
        if node is self :
            return False
        if  (self.net.maxparent and len (self.parents) >= self.net.maxparent) :
            return False
        if self.is_transitive_parent (node) :
            return False
        return True
    # end def may_append_parent

    def __hash__ (self) :
        return self.idx
    # end def __hash__

    def __repr__ (self) :
        if self.parents :
            self.rank = max (self.rank, max (p.rank for p in self.parents) + 1)
        r = [ "%sNode: %s children: %s parents: %s"
            % ( '-' * self.rank
              , self.idx
              , tuple (c.idx for c in self.children)
              , tuple (p.idx for p in self.parents)
              )
            ]
        if self.verbose :
            indent = ' ' * self.rank
            dn = str (self.dnode)
            dn = '\n'.join (indent + d for d in dn.split ('\n'))
            r.append (dn)
        #for c in self.children :
        #    r.append (str (c))
        return '\n'.join (r)
    # end def __repr__
    __str__ = __repr__

# end class BNode

class DNode (object) :
    """ Binary decision tree node
        A Tree either has two children.
        A child may be a leaf not (which contains a single probability).
        or another DNode.
    """

    def __init__ (self, bnode, parent, cidx = None) :
        self.bnode    = bnode
        self.debug    = self.bnode.debug
        self.idx      = bnode.idx
        self.children = []
        self.parent   = parent
        if cidx is None :
            self.genes = parent.genes
        else :
            self.genes = parent.gsplit [cidx]
        self.lvl      = self.parent.lvl + 1
        self.gsplit   = [[], []]
        self.n        = 0
        for g in self.genes :
            self.gsplit [g [self.idx]].append (g)
            self.n += 1
    # end def __init__

    def __repr__ (self) :
        indent = '  ' * self.lvl
        r = []
        r.append ("%s%d" % (indent, self.idx))
        for c in self.children :
            r.append (str (c))
        return '\n'.join (r)
    # end def __repr__
    __str__ = __repr__

    def feasible (self) :
        assert self.idx == self.bnode.idx
        assert self.idx != self.parent.idx
        f = self.children [0].bnode.may_append_parent (self.bnode)
        self.debug \
            ( "feasible %s: split %s on %s"
            % (f, self.children [0].bnode.idx, self.idx)
            )
        return f
    # end def feasible

# end class DNode

class DLeaf (object) :
    """ Binary decision tree leaf
    """

    def __init__ (self, bnode, parent, cidx, genes) :
        self.bnode      = bnode
        self.idx        = bnode.idx
        self.cidx       = cidx
        self.parent     = parent
        self.lvl        = parent.lvl + 1
        self.genes      = genes
        self.candidates = {}
        self.by_gain    = None
        self.debug      = self.bnode.debug
        self.min_n      = self.bnode.net.min_n
        self.n  = 0
        self.n1 = 0
        for g in self.genes :
            if g [self.idx] :
                self.n1 += 1
            self.n += 1
        if self.n == 0 :
            self.p = 1.0
        else :
            self.p = 1.0 * self.n1 / self.n
        lgamma = self.bnode.net.lgamma
        self.score = 0.0
        self.score += lgamma (1 + self.n1)
        self.score += lgamma (1 + self.n - self.n1)
        self.score -= lgamma (2 + self.n)
        if not isinstance (self.parent, BNode) :
            self.parent.children.append (self)
            assert len (self.parent.children) <= 2
    # end def __init__

    def __repr__ (self) :
        indent = '  ' * self.lvl
        return ("%s%d: %1.4f" % (indent, self.idx, self.p))
    # end def __repr__
    __str__ = __repr__

    def candidate_iter (self) :
        """ Iterate over candidates in sorted (best first) order
            We take care to check if the candidate is still existing
        """
        if self.by_gain is None :
            self.by_gain = list \
                (sorted
                  ( self.candidates.keys ()
                  , key = lambda x : -self.candidates [x].gain
                  )
                )
        for idx in self.by_gain :
            if idx in self.candidates :
                yield (self.candidates [idx])
    # end def candidate_iter

    def del_candidate (self, cand) :
        del self.candidates [cand.idx]
    # end def del_candidate

    def try_add_candidate (self, bnode) :
        """ Try split on bnode
        """
        if self.n < self.min_n :
            return 0.0
        idx = bnode.idx
        assert idx not in self.candidates
        assert idx != self.idx
        p = self.parent
        # Do we already have a split on that idx?
        while not isinstance (p, BNode) :
            if p.idx == idx :
                return
            p = p.parent
        cidx = self.cidx
        if isinstance (self.parent, BNode) :
            cidx = None
        n  = DNode (bnode, self.parent, cidx)
        c1 = self.__class__ (self.bnode, n, 0, genes = n.gsplit [0])
        c2 = self.__class__ (self.bnode, n, 1, genes = n.gsplit [1])
        n.gain = c1.score + c2.score - self.score - self.bnode.net.cutoff
        if n.gain > 0 :
            self.candidates [n.idx] = n
        return n.gain
    # end def try_add_candidate

    def split (self, idx) :
        cand = self.candidates [idx]
        if isinstance (self.parent, BNode) :
            assert self.cidx == 0
            self.parent.dnode = cand
        else :
            self.parent.children [self.cidx] = cand
        self.debug ("split %s on %s" % (self.idx, cand.bnode.idx))
        self.debug ("split (dnode):", cand.idx, self.parent.idx, end = ' ')
        self.debug ("split (bnode):", cand.bnode.idx, self.bnode.idx, end = ' ')
        self.debug ("cbnode:", cand.children [0].idx)
        self.bnode.append_parent (cand.bnode)
        for node in self.bnode.net.nodes :
            if self.bnode.may_append_parent (node) :
                for l in cand.children :
                    g = l.try_add_candidate (node)
    # end def split

# end class DLeaf

class HBOA (PMBGA) :
    """ hierarchical Bayesian Optimization Algorithm
    """

    def build_model (self, p_pop) :
        self.__super.build_model (p_pop)
        self.net = Bayesian_Network (self, self.genes, lgamma = self.lgamma)
    # end def build_model

    def clear_cache (self) :
        pass
    # end def clear_cache

    def lgamma (self, x) :
        """ lgamma cache, about 3 times faster than calling lgamma
        """
        return self.l_gamma [x]
    # end def lgamma

    def sample_individual (self, p, pop) :
        d = self.net.sample_model ()
        assert len (d) == len (self)
        for k in sorted (d) :
            self.set_allele (p, pop, k, d [k])
    # end def sample_individual

    def post_init (self) :
        self.__super.post_init ()
        self.do_debug = 0
        self.verbose  = 0
        self.clear_cache ()
        self.l_gamma = [0]
        for k in range (self.pop_size + 2) :
            self.l_gamma.append (lgamma (k + 1))
    # end def post_init

    def print_model (self) :
        for r in self.net.nodes :
            print (r, file = self.file)
    # end def print_model

# end class HBOA
