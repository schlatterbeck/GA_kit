+++++++++++++++++++++++++++
Advanced Genetic Algorithms
+++++++++++++++++++++++++++

This package implements some of the advanced algorithms on top of my
python wrapper ``pgapy`` of the Parallel Genetic Algorithm Package
PGAPack. In the following we're using GA to mean genetic algorithm.

Currently we cover two variants of probabilistic model building GAs,
also called Estimation of Distribution Algorithm (EDA) by other authors.
The variants implemented here all operate on binary genes.

You can implement your own solver for an optimization problem on top of
this implementation by studying the documentation of pgapy_ and the user
guide of PGAPack_ (which comes bundled with an installation of pgapy_).
In addition looking into the implementation of the deceptive function it
becomes clear how to code your own problem.

Deceptive Functions
===================

For testing these algorithms, deceptive functions have been proposed
[DG92]_, these work on *k* bit strings and all, except one value lead
the genetic algorithm in the wrong direction. A concrete example are *k*
bit trap functions. The function evaluates to the number of 1-bits minus
1, *except* for the single case with all zeros for which the evaluation
is *k*. So for all values the algorithm gains when adding more 1-bits
but the case with all 0-bits is the optimum. Several of these trap
functions can be concatenated into a more complex problem, the
evaluation being the sum of all individual trap functions.

An implementation of order *k* trap functions is given in
``deceptive.py``, the number of bits and the number of concatenated trap
functions can be specified on the command-line. By default, all trap
functions are concatenated in order, making the problem solveable by a
simple GA. With the ``--shuffle`` option, the bits belonging to
individual trap functions are shuffled across the whole gene. Due to the
disruptive effect of crossover in the simple GA, with shuffled genes,
the simple GA can no longer solve these problems.

The newer algorithms build a probabilistic model during optimization.
This model can detect which genes have a high correlation and should be
treated as a unit. This makes even shuffled deceptive problems solveable
by these algorithms.

Probabilistic Model Building Genetic Algorithms (PMBGA)
=======================================================

This implements some of the probabilistic model building genetic
algorithms from the literature. We also implement the typical test
functions the so-called deceptive functions. As a basis it needs PGApy_,
my wrapper around the Parallel Genetic Algorithm Library, PGApack_.
For the latter I've recently implemented restricted tournament
replacement [Har95]_ as one of the replacement strategies which is better at
preserving genetic variability for longer that the standard replacements
schemes.

PMBGAs typically replace the mutation and crossover operators of
traditional genetic algorithms [Hol75]_, [Gol89]_ with a statistical
model built during the search. Instead of crossover and mutation the
model is sampled to create the next generation.

Extended Compact Genetic Algorithm (ECGA)
-----------------------------------------

The ECGA, originally proposed by Harik [Har99]_ replaces the mutation and
crossover of a genetic algorithm by building and later sampling
a statistical model of the population in each generation. The
model uses a minimum description length (MDL) model. A better
description with an example is found in [HLS06]_.
Note that we're using the restricted tournament replacement,
originally called restricted tournament selection by Harik
[Har94]_, [Har95]_ and later named restricted tournament replacement by
Pelikan [Pel05]_ for population replacement which has recently been
implemented in the genetic algorithm library PGApack.
Note that currently we support only bitstring (bool) genes.
Since we're using model building not the standard GA operations,
we disable mutation by setting the mutation probability to 0.

Hierarchical Bayesian Optimization Algorithm (hBOA)
---------------------------------------------------

Bayesian optimization, first introduced as the Bayesian Optimization
Algorithm (BOA) [PGC99]_ uses a Bayesian network to model the
dependencies among the individual genes. It was later extended to
include local structures in the Bayesian network which allows more
granular model building and named hBOA [Pel05]_.  The hBOA, originally
proposed by Pelikan and Goldberg [PG00]_, [PG01]_ and described in Pelikan's
dissertation [Pel02]_ which was developed in a book on Bayesian
optimization algorithms [Pel05]_ builds a Bayesian network to model the
dependencies of the variables in the genetic algorithm. Like other
PMBGAs it samples the resulting model to generate individuals as
candidates for the next generation. It uses restricted tournament
replacement (RTR) [Har95]_ as the replacement strategy. For implementing
it using PGApy, we use the RTR strategy from PGApack which was recently
introduced. When building the model it uses a Bayesian dirichlet metric
and a term that penalizes models which are too complex ([Pel05]_ p.113).
In the code, this term is named *cutoff* because it stops the addition
of further edges to the Bayesian network. It was later discovered that
when using noisy selection strategies like tournament selection, the 
algorithm would add spurious dependencies [LLPG10]_, [LLPG11]_,
especially during the
final phase of each network-building stage. This could be countered by
adding a penalty factor to the *cutoff* parameter. The best value
according to the tests was the tournament size, e.g. 2 for binary
tournament, this is the reason the factor was called *s-penalty*. To
further limit the spurious dependencies generated the code has an
additional parameter ``min_split`` that will not split the graph further
when the number of samples affected is below the threshold given by the
``min_split`` parameter.

.. [DG92] Kalyanmoy Deb and David E. Goldberg. Analyzing deception in
   trap functions. In L. Darrell Whitley, editor, Foundation of Genetic
   Algorithms (FOGA) 2, pages 93–108.  Elsevier, 1992.
.. [Gol89] David E. Goldberg. Genetic Algorithms in Search, Optimization
   & Machine Learning. Addison Wesley, October 1989.
.. [Har94] Georges R. Harik. Finding multiple solutions in problems of
   bounded difficulty. IlliGAL Report 94002, Illinois Genetic
   Algorithm Lab, May 1994.
.. [Har95] Georges R. Harik. Finding multimodal solutions using
   restricted tournament selection. In Larry J. Eshelman, editor,
   Proceedings of the International Conference on Genetic Algorithms
   (ICGA), pages 24–31. Morgan Kaufmann, July 1995.
.. [Har99] George R. Harik. Linkage learning via probabilistic modeling
   in the ECGA. IlliGAL Report 99010, Illinois Genetic Algorithm Lab,
   January 1999
.. [HLS06] Georges R. Harik, Fernando G. Lobo, and Kumara Sastry.
   Linkage learning via probabilistic modeling in the extended compact
   genetic algorithm (ECGA). In Martin Pelikan, Kumara Sastry, and
   Erick Cantú-Paz, editors, Scalable Optimization via Probabilistic
   Modelling: From Algorithms to Applications, volume 33 of Studies in
   Computational Intelligence, pages 39–61. Springer, 2006.
.. [Hol75] John Holland. Adaptation in Natural and Artificial Systems.
   University of Michigan Press, Ann Arbor, Michigan, 1975.
.. [LLPG10] Claudio F. Lima, Fernando G. Lobo, Martin Pelikan, and David
   E. Goldberg. Model accuracy in the Bayesian optimization algorithm.
   IlliGAL Report 2010002, Illinois Genetic Algorithm Lab, March 2010.
.. [LLPG11] Claudio F. Lima, Fernando G. Lobo, Martin Pelikan, and David
   E. Goldberg. Model accuracy in the Bayesian optimization algorithm.
   Soft Computing, 15(7):1351–1371, July 2011.
.. [Pel02] Martin Pelikan. Bayesian Optimization Algorithm: From Single
   Level to Hierarchy. Dissertation, University of Michigan, 2002.
   IlliGAL Report No. 2002023.
.. [Pel05] Martin Pelikan. Hierarchical Bayesian Optimization
    Algorithm: Toward a New Generation of Evolutionary
    Algorithms, volume 170 of Studies in Fuzziness and Soft
    Computing. Springer, 2005.
.. [PG00] Martin Pelikan and David E. Goldberg. Hierarchical problem
   solving by the Bayesian optimization algorithm. In L. Darrell Whitley,
   David E. Goldberg, Erick Cantú-Paz, Lee Spector, Ian C. Parmee, and
   Hans-Georg Beyer, editors, Genetic and Evolutionary Computation
   GECCO 2000, pages 267–274, Las Vegas, Nevada, July 2000. Morgan
   Kaufmann.
.. [PG01] Martin Pelikan and David E. Goldberg. Escaping hierarchical
   traps with competent genetic algorithms. In Lee Spector, Erik D.
   Goodman, Annie Wu, William B. Langdon, and Hans Michael Voigt,
   editors, Genetic and Evolutionary Computation Conference
   (GECCO-2001), pages 511–518, Seattle, WA, July 2001. Morgan Kaufmann.
.. [PGC99] Martin Pelikan, David E. Goldberg, and Erick Cantú-Paz. BOA:
   The Bayesian optimization algorithm. In Wolfgang Banzhaf, Jason M.
   Daida, A. E. Eiben, Max H. Garzon, Vasant G. Honavar, Mark J.
   Jakiela, and Robert E. Smith, editors, Genetic and Evolutionary
   Computation GECCO 1999, page 525–532, Orlando, Florida, July 1999.
   Morgan Kaufmann.

.. _PGApy: https://github.com/schlatterbeck/pgapy
.. _PGApack: https://github.com/schlatterbeck/pgapack
