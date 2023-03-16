This repository contains the Matlab code used for the paper ["A posteriori validation of generalized polynomial chaos expansions"](https://arxiv.org/abs/2203.02404) by M. Breden.

In order to reproduce the computer-assisted parts of the proofs discussed in the paper, download the whole repository and run
* script_BasicExample.m for Proposition 2.16, 
* script_Lorenz.m  for Theorem 3.5, and to reproduce the comparisons presented in Table 2 and Table 3 (section 3.4),
* script_SwiftHohenberg.m for Theorem 4.1,
* script_SwiftHohenberg_2para.m for Theorem 4.3.

You can also use the other scripts (script_LorenzExplore.m, script_SwiftHohenberg_continuation.m and script_SwiftHohenberg_Explore.m), to study different solutions and parameter values than those discussed in the paper.

The code does not require the Intlab toolbox (http://www.ti3.tu-harburg.de/intlab/) to run, but without it the results are not guaranteed to be immune to rounding errors.
