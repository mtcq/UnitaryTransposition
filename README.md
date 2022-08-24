## Code to accompany: [Experimental superposition of time directions](https://arxiv.org/abs/xxxx.xxxx)

#### Teodor Strömberg, Peter Schiansky, Marco Túlio Quintino, Michael Antesberger, Lee Rozema, Iris Agresti, Caslav Brukner, and Philip Walther


This is a repository for the article [Experimental superposition of time directions](https://arxiv.org/abs/xxxx.xxxx)

 This code requires:
- [cvx](http://cvxr.com/) - a free MATLAB toolbox for rapid prototyping of optimization problems.
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory.

This repository consists of:

- [MaxSuccessProb_Utrans.m](https://github.com/mtcq/UnitaryTransposition/blob/main/MaxSuccessProb_Utrans.m):
Script finds the maximal success probabilities for Parallel, Sequential, and general strategies with standard SDP methods.

- [MaxSuccessProb_Utrans_symbolic.m](https://github.com/mtcq/UnitaryTransposition/blob/main/MaxSuccessProb_Utrans_symbolic.m):
Script finds and verifies the upper bounds for all strategy classes in a rigorous computer assisted proof way.

- [DefineSets_Utrans.m](https://github.com/mtcq/UnitaryTransposition/blob/main/DefineSets_Utrans.m):
Script defines the sets considered in the discrimination task.


- [DefineSets_Utrans_symbolic.m](https://github.com/mtcq/UnitaryTransposition/blob/main/DefineSets_Utrans_symbolic.m):
Script defines the sets considered in the discrimination task as symbolic variables

- [auxiliary_functions](https://github.com/mtcq/UnitaryTransposition/tree/main/auxiliary_functions):
Folder with auxiliary functions.
