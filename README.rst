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

Extended Compacat Genetic Algorithm (ECGA)
------------------------------------------

The ECGA, originally proposed in [Har99]_ replaces the mutation and
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
Algorithm (BOA) [PGC99]_ uses a bayesian network to model the
dependencies among the individual genes. It was later extended to
include local structures in the bayesian network which allows more
granular model building and named hBOA [Pel05]_.  The hBOA, originally
proposed by Pelikan and Goldberg [PG00]_, [PG01]_ and described in Pelikan's
dissertation [Pel02]_ which was developed in a book on bayesian
optimizatino algorithms [Pel05]_ builds a bayesian network to model the
dependencies of the variables in the genetic algorithm. Like other
PMBGAs it samples the resulting model to generate individuals as
candidates for the next generation. It uses restricted tournament
replacement (RTR) [Har95]_ as the replacement strategy. For implementing
it using PGApy, we use the RTR strategy from PGApack which was recently
introduced. When building the model it uses a bayesian dirichlet metric
and a term that penalizes models which are too complex ([Pel05]_ p.113).
In the code, this term is named *cutoff* because it stops the addition
of further edges to the bayesian network. It was later discovered that
when using noisy selection strategies like tournament selection, the 
algorithm would add spurious dependencies [LLPG10]_, especially during the
final phase of each network-building phase. This could be countered by
adding a penalty factor to the *cutoff* parameter. The best value
according to the tests was the tournament size, e.g. 2 for binary
tournament, this is the reason the factor was called *s-penalty*. To
further limit the spurious dependencies generated the code has an
additional parameter ``min_split`` that will not split the graph further
when the number of samples affected is below the threshold given by the
``min_split`` parameter.

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
   E. Goldberg. Model accuracy in the bayesian optimization algorithm.
   IlliGAL Report 2010002, Illinois Genetic Algorithm Lab, March 2010.
.. [Pel02] Martin Pelikan. Bayesian Optimization Algorithm: From Single
   Level to Hierarchy. Dissertation, University of Michigan, 2002.
   IlliGAL Report No. 2002023.
.. [Pel05] Martin Pelikan. Hierarchical Bayesian Optimization
    Algorithm: Toward a New Generation of Evolutionary
    Algorithms, volume 170 of Studies in Fuzziness and Soft
    Computing. Springer, 2005.
.. [PG00] Martin Pelikan and David E. Goldberg. Hierarchical problem
   solving by the bayesian optimization algorithm. In L. DarrellWhitley,
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
   The bayesian optimization algorithm. In Wolfgang Banzhaf, Jason M.
   Daida, A. E. Eiben, Max H. Garzon, Vasant G. Honavar, Mark J.
   Jakiela, and Robert E. Smith, editors, Genetic and Evolutionary
   Computation GECCO 1999, page 525–532, Orlando, Florida, July 1999.
   Morgan Kaufmann.

.. _PGApy: https://github.com/schlatterbeck/pgapy
.. _PGApack: https://github.com/schlatterbeck/pgapack
