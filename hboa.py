#!/usr/bin/python3

from sys  import stderr
from math import log, lgamma
from sga  import PMBGA, log2

class BNode (object) :

    def __init__ (self, hboa, idx) :
        self.hboa     = hboa
        self.debug    = hboa.debug
        self.genes    = hboa.genes
        self.idx      = idx
        self.parents  = {}
        self.children = {}
        self.lvl      = 0
        self.dnode    = DLeaf (self, self, 0, self.genes)
        self.rank     = 0
    # end def __init__

    def add_initial_candidates (self) :
        for idx in self.hboa.nodes :
            if idx == self.idx :
                continue
            self.dnode.try_add_candidate (self.hboa.nodes [idx])
    # end def add_initial_candidates

    def append_parent (self, node) :
        self.debug ("append_parent: %d %d" % (self.idx, node.idx))
        assert node not in self.children
        assert self not in node.parents
        self.parents  [node] = 1
        node.children [self] = 1
        if self in self.hboa.roots :
            del self.hboa.roots [self]
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
        for p in self.parents :
            if node is p :
                return True
            if p.is_transitive_parent (node) :
                return True
        return False
    # end def is_transitive_parent

    def is_transitive_child (self, node) :
        for c in self.children :
            if node is c :
                return True
            if c.is_transitive_child (node) :
                return True
        return False
    # end def is_transitive_child

    def may_append_parent (self, node) :
        if node is self :
            return False
        if self.is_transitive_parent (node) :
            return False
        if self.is_transitive_child (node) :
            return False
        return True
    # end def may_append_parent

    def __hash__ (self) :
        return self.idx
    # end def __hash__

    def __repr__ (self) :
        if self.parents :
            self.rank = max (self.rank, max (p.rank for p in self.parents) + 1)
        indent = '  ' * self.rank
        r = [ "%sNode: %s children: %s parents: %s"
            % ( '--' * self.rank
              , self.idx
              , tuple (c.idx for c in self.children)
              , tuple (p.idx for p in self.parents)
              )
            ]
        dn = str (self.dnode)
        dn = '\n'.join (indent + d for d in dn.split ('\n'))
        r.append (dn)
        for c in self.children :
            r.append (str (c))
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

    def __init__ (self, bnode, parent) :
        self.bnode    = bnode
        self.debug    = self.bnode.debug
        self.idx      = bnode.idx
        self.children = []
        self.parent   = parent
        self.genes    = parent.genes
        self.lvl      = self.parent.lvl + 1
        self.g0       = []
        self.g1       = []
        self.n1       = 0
        self.n        = 0
        for g in self.genes :
            if g [self.idx] :
                self.n1 += 1
                self.g1.append (g)
            else :
                self.g0.append (g)
            self.n += 1
        self.p = 1.0 * self.n1 / self.n
    # end def __init__

    def __repr__ (self) :
        indent = '  ' * self.lvl
        r = []
        r.append ("%s%d: %1.4f" % (indent, self.idx, self.p))
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
        self.debug      = self.bnode.debug
        self.n  = 0
        self.n1 = 0
        for g in self.genes :
            if g [self.idx] :
                self.n1 += 1
            self.n += 1
        self.p = 1.0 * self.n1 / self.n
        lgamma = self.bnode.hboa.lgamma
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

    def try_add_candidate (self, bnode) :
        """ Try split on bnode
        """
        idx = bnode.idx
        assert idx not in self.candidates
        assert idx != self.idx
        p = self.parent
        # Do we already have a split on that idx?
        while not isinstance (p, BNode) :
            if p.idx == idx :
                return
            p = p.parent
        n  = DNode (bnode, self.parent)
        c1 = self.__class__ (self.bnode, n, 0, genes = n.g0)
        c2 = self.__class__ (self.bnode, n, 1, genes = n.g1)
        n.gain = c1.score + c2.score - self.score - self.bnode.hboa.l2n2
        if n.gain > 0 :
            self.candidates [n.idx] = n
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
        for x in self.bnode.hboa.nodes :
            node = self.bnode.hboa.nodes [x]
            if self.bnode.may_append_parent (node) :
                for l in cand.children :
                    l.try_add_candidate (node)
    # end def split

# end class DLeaf

class HBOA (PMBGA) :
    """ hierarchical Bayesian Optimization Algorithm
    """

    def build_model (self, p_pop) :
        self.__super.build_model (p_pop)
        self.nodes = {}
        self.roots = {}
        for bitidx in range (len (self)) :
            node = BNode (self, bitidx)
            self.nodes [node.idx] = node
            self.roots [node] = 1
        # The initial candidates already need the full list of nodes
        for bitidx in self.nodes :
            self.nodes [bitidx].add_initial_candidates ()
        while 1 :
            maxgain = -1
            leave   = None
            cidx    = -1
            for idx in self.nodes :
                node = self.nodes [idx]
                for l in node.leaves () :
                    keys = list (l.candidates.keys ())
                    for idx in keys :
                        c = l.candidates [idx]
                        if c.feasible () :
                            if c.gain > maxgain :
                                maxgain = c.gain
                                leave   = l
                                cidx    = idx
                        else :
                            del l.candidates [idx]
            if maxgain <= 0 :
                break
            leave.split (cidx)
    # end def build_model

    def clear_cache (self) :
        pass
    # end def clear_cache

    def debug (self, *args, **kw) :
        if self.do_debug :
            print (*args, **kw, file = stderr)
    # end def debug

    def lgamma (self, x) :
        """ lgamma cache, about 3 times faster than calling lgamma
        """
        return self.l_gamma [x]
    # end def lgamma

    def sample_individual (self, p, pop) :
        d = {}
        for r in self.roots :
            c = r.dnode
            assert isinstance (c, DLeaf)
    # end def sample_model

    def post_init (self) :
        self.__super.post_init ()
        self.do_debug = False
        self.clear_cache ()
        self.l2n2 = log2 (len (self)) / 2.0
        self.l_gamma = [0]
        for k in range (self.pop_size + 2) :
            self.l_gamma.append (lgamma (k + 1))
    # end def post_init

    def print_model (self) :
        for r in self.roots :
            print (r, file = self.file)
    # end def print_model

# end class HBOA
