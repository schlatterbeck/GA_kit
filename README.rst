Probabilistic Model Building Genetic Algorithms (PMBGA)
=======================================================

This implements some of the probabilistic model building genetic
algorithms from the literature. We also implement the typical test
functions the so-called deceptive functions. As a basis it needs PGApy_,
my wrapper around the Parallel Genetic Algorithm Library, PGApack_.
For the latter I've recently implemented restricted tournament
replacement [5]_ as one of the replacement strategies which is better at
preserving genetic variability for longer that the standard replacements
schemes.

Extended Compacat Genetic Algorithm (ECGA)
------------------------------------------

The ECGA, originally proposed in [1]_ replaces the mutation and
crossover of a genetic algorithm by building and later sampling
a statistical model of the population in each generation. The
model uses a minimum description length (MDL) model. A better
description with an example is found in [2]_.
Note that we're using the restricted tournament replacement,
originally called restricted tournament selection by Harik
[4]_, [5]_ and later named restricted tournament replacement by
Pelikan [7]_ for population replacement which has recently been
implemented in the genetic algorithm library PGApack.
Note that currently we support only bitstring (bool) genes.
Since we're using model building not the standard GA operations,
we disable mutation by setting the mutation probability to 0.

.. [1] George R. Harik. Linkage learning via probabilistic modeling in
    the ECGA. IlliGAL Report 99010, Illinois Genetic Algorithm
    Lab, January 1999
.. [2] Georges R. Harik, Fernando G. Lobo, and Kumara Sastry.
    Linkage learning via probabilistic modeling in the extended
    compact genetic algorithm (ECGA). In Pelikan et al. [3]_,
    pages 39-61.
.. [3] Martin Pelikan, Kumara Sastry, and Erick Cantu-Paz, editors.
    Scalable Optimization via Probabilistic Modelling: From
    Algorithms to Applications, volume 33 of Studies in
    Computational Intelligence. Springer, 2006.
.. [4] Georges R. Harik. Finding multiple solutions in problems of
    bounded difficulty. IlliGAL Report 94002, Illinois Genetic
    Algorithm Lab, May 1994.
.. [5] Georges R. Harik. Finding multimodal solutions using
    restricted tournament selection. In Eshelman [6]_, pages
    24-31.
.. [6] Larry J. Eshelman, editor. Proceedings of the International
    Conference on Genetic Algorithms (ICGA). Morgan Kaufmann
    July 1995.
.. [7] Martin Pelikan. Hierarchical Bayesian Optimization
    Algorithm: Toward a New Generation of Evolutionary
    Algorithms, volume 170 of Studies in Fuzziness and Soft
    Computing. Springer, 2005.

.. _PGApy: https://github.com/schlatterbeck/pgapy
.. _PGApack: https://github.com/schlatterbeck/pgapack
