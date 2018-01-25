=======
Cutnorm
=======

Approximation via Gaussian Rounding and Optimization with Orthogonality Constraints
-----------------------------------------------------------------------------------

This package computes the approximations to the cutnorm using some of the techniques detailed by Alon and Noar [ALON2004]_ and a fast optimization algorithm by Wen and Yin [WEN2013]_.

Read the documentation_.

.. _documentation: https://pingkoc.github.io/cutnorm/cutnorm.html

Installation
------------

Use pip_ to install the package.
Install from terminal as follows::

  $ pip install cutnorm

.. _pip: http://www.pip-installer.org/en/latest/

Example Usage
-------------

Below is an example of using the cutnorm package and tools. Given two graphs A and B, we wish to compute a norm for the difference matrix (A - B) between the two graphs. An obvious example to represent the advantage of using a cutnorm over l1 norm is to consider A and B as `Erdos-Renyi random graphs`_. Under a fixed vertex set, an Erdos-Renyi random graph is one where a fixed probability determines the presence of an edge.

.. _`Erdos-Renyi random graphs`: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model

Given two Erdos-Renyi random graphs with fix n and p=0.5, the edit distance (l1 norm) of the difference (after normalization) is 1/2 with large probability. However, these two graphs have the same global structure. The edit distance fails as a notion of 'distance' between the two graphs in the perspective of global structural similarity as discussed by Lovasz [LOVASZ2009]_. The cutnorm is a measure of distance that reflects global structural similarity. In fact, the cutnorm of the difference for this example approaches 0 as n grows.

.. code:: python

    import numpy as np
    from cutnorm import cutnorm, tools

    # Generate Erdos Renyi Random Graph
    n = 100
    p = 0.5
    erdos_renyi_a = tools.sbm.erdos_renyi(n, p)
    erdos_renyi_b = tools.sbm.erdos_renyi(n, p)

    # Compute l1 norm
    normalized_diff = (erdos_renyi_a - erdos_renyi_b) / n**2
    l1 = np.linalg.norm(normalized_diff.flatten(), ord=1)

    # Compute cutnorm
    cutn_round, cutn_sdp, info = cutnorm(erdos_renyi_a, erdos_renyi_b)

    print("l1 norm: ", l1)                 # prints l1 norm value near ~0.5
    print("cutnorm rounded: ", cutn_round) # prints cutnorm rounded solution near ~0
    print("cutnorm sdp: ", cutn_sdp)       # prints cutnorm sdp solution near ~0

----

.. [ALON2004] Noga Alon and Assaf Naor. 2004. Approximating the cut-norm via Grothendieck's inequality. In Proceedings of the thirty-sixth annual ACM symposium on Theory of computing (STOC '04). ACM, New York, NY, USA, 72-80. DOI: http://dx.doi.org/10.1145/1007352.1007371
.. [WEN2013] Zaiwen Wen and Wotao Yin. 2013. A feasible method for optimization with orthogonality constraints. Math. Program. 142, 1-2 (December 2013), 397-434. DOI: https://doi.org/10.1007/s10107-012-0584-1
.. [LOVASZ2009] Lovasz, L. 2009. Very large graphs. ArXiv:0902.0132 [Math]. Retrieved from http://arxiv.org/abs/0902.0132
