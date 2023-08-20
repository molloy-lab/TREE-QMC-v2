TREE-QMC
========

TREE-QMC is a quartet-based method for estimating species trees from gene trees. The input/output are the same as the popular method ASTRAL. To learn more about TREE-QMC, check out the related [paper](http:doi.org/10.1101/gr.277629.122) in *Genome Research* and [abstract](https://link.springer.com/book/10.1007/978-3-031-29119-7) in the *RECOMB 2023 proceedings*.

Note that this repository corresponds to the latest version of TREE-QMC (version 2). A preprint will be available in the near future with details. Note that the original code (version 1) is available at [this Github repository](https://github.com/molloy-lab/TREE-QMC-v1).

Acknowledgements
-----------------
TREE-QMC builds upon the algorithmic framework work introduced by Sagi Snir and Satish Rao; see [this *IEEE/ACM TCBB* paper](http:doi.org/10.1109/TCBB.2008.133) as well as [this *Systematic Biology* paper](http:doi.org/10.1093/sysbio/syu087).

Usage
-----

To build, TREE-QMC use commands:
```
git clone https://github.com/molloy-lab/TREE-QMC.git
cd TREE-QMC/MQLib
make
cd ..
g++ -std=c++11 -O2 -I MQLib/include -o TREE-QMC src/*.cpp MQLib/bin/MQLib.a
```

To run TREE-QMC, use command:
```
./TREE-QMC -i <input file> -o <output file name>
```
See [this tutorial](example/tutorial.md) for more information. Also, check out the TREE-QMC usage options with this command:
```
./TREE-QMC -h
```

The output should be
```
TREE-QMC version 2.0.0
COMMAND: ./TREE-QMC -h 
=================================== TREE-QMC ===================================
This is version 2.0.0 of TREe Embedded Quartet Max Cut (TREE-QMC).

USAGE:
./TREE-QMC (-i|--input) <input file> [(-o|--output) <output file>]
           [(--polyseed) <integer>] [(--maxcutseed) <integer>]
           [(-n|--normalize) <normalization scheme>]
           [(-x|--execution) <execution mode>]
           [(-v|--verbose) <verbose mode>] [-h|--help]

OPTIONS:
[-h|--help]
        Prints this help message.
(-i|--input) <input file>
        Name of file containing gene trees in newick format (required)
        IMPORTANT: current implementation of TREE-QMC requires that the input
        gene trees are unrooted and binary. Thus, TREE-QMC suppresses roots
        and randomly refines polytomies during a preprocessing phase; the
        resulting trees are written to "<input file>.refined".
[(-o|--output) <output file>]
        Name of file for writing output species tree (default: stdout)
[(--polyseed) <integer>]
        Seeds random number generator with <integer> prior to arbitrarily
        resolving polytomies. If <integer> is set to -1, system time is used;
        otherwise, <integer> should be positive (default: 12345).
[(--maxcutseed) <integer>]
        Seeds random number generator with <integer> prior to calling the max
        cut heuristic but after the preprocessing phase. If <integer> is set to
        -1, system time is used; otherwise, <integer> should be positive
        (default: 1).
[(-n|--normalize) <normalization scheme>]
        Initially, each quartet is weighted by the number of input gene
        trees that induce it. At each step in the divide phase of wQMC and
        TREE-QMC, the input quartets are modified with artificial taxa. We
        introduce two normalization schemes for artificial taxa and find
        that they improve empirical performance of TREE-QMC in a simulation
        study. The best scheme is run by default. See paper for details.
        -n 0: none
        -n 1: uniform
        -n 2: non-uniform (default)
[(-x|--execution) <execution mode>]
        TREE-QMC uses an efficient algorithm that operates directly on the
        input gene trees by default. The naive algorithm, which operates on a
        set of quartets weighted based on the input gene trees, is also
        implemented for testing purposes.
        -x 0: run efficient algorithm (default)
        -x 1: run naive algorithm
[(-v|--verbose) <verbose mode>]
        -v 0: write no subproblem information (default)
        -v 1: write CSV with subproblem information (subproblem ID, parent
              problem ID, depth of recursion, total # of taxa, # of artifical
              taxa, species names)
        -v 2: write CSV with subproblem information (info from v1 plus # of
              of elements in f, # of pruned elements in f, # of zeroes in f)
[--shared <use shared taxon data structure to normalize quartet weights>]
        Do NOT use unless there are no missing data!!!

Contact: Post issue to Github (https://github.com/molloy-lab/TREE-QMC/)
        or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)

If you use TREE-QMC in your work, please cite:
  Han and Molloy, 2023, "Improving quartet graph construction for scalable
  and accurate species tree estimation from gene trees," Genome Research,
  http:doi.org/10.1101/gr.277629.122.
================================================================================
```
