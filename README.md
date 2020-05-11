# GKC: Graph Kernel Collection
Public Release of parts of CMU-SPEED's implementations of GAP benchmark suite.
The 6 graph algorithms in the GAP Benchmark Suite are represented here (BFS, 
Betweeness Centrality, Connected Components, Pagerank, SSSP, and Triangle Counting).

## How to run
The top-level directory for each algorithm contains a base file with the
benchmark driving code, as well as a static library (\*.a) containing compiled
code for each algorithm implementation. These libraries were compiled for
Skylake-X CPUs (on an Intel Xeon Platinum 8153) with flags for AVX2 and AVX-512
enabled.

The verification and non-verification binaries can be compiled from each
algorithm directory by running ```make```. Most algorithms will compile two 
executables - one with and one without verification. 
See the compilation section for more details on compilation setup.

### Compilation
See the makefiles in each algorithm's directory for compile-time flags used.
The main driver file (typically main.cpp, or [name of algorithm].cpp) can be
modified and recompiled with the static library provided for the algorithm.

The static library file (.a) in each algorithm directory has been compiled with
icpc (ICC) version 19.1 20200306. For best results, the remainder of the code
should be compiled with the same compiler.

To use icpc (ICC) 19.1 20200306 on devcloud, you should add the following lines
to your bash\_profile:
```
# Enable Intel tools
export INTEL_LICENSE_FILE=/usr/local/licenseserver/psxe.lic
export PATH=/glob/intel-python/python3/bin/:/glob/intel-python/python2/bin/:${PATH}
source /glob/development-tools/parallel-studio/bin/compilervars.sh intel64
export PATH=$PATH:/bin
if [ -d /opt/intel/inteloneapi ]; then source /opt/intel/inteloneapi/setvars.sh > /dev/null 2>&1; fi
```

### PBS Scripts
Each algorithm directory has a sub-folder called pbs/. This directory contains
sample portable batch scripts (PBS) for running the algorithm with the 5 GAP
graphs. By default the executable used is in the pbs/ directory.
Outputs from the pbs scripts go into .dat files to the pbs/outputs/ subdirectory.

The PBS scripts contain the OMP runtime parameters we used for testing our
implementations on 32 and 64 threads on a Platinum 8153 system, using
interleaved memory allocation via numactl and thread pinning to control for SMT
and non-SMT operation.

If you intend to use the PBS scripts, you will need to update the base
directory, data directory, and possibly sources directory in each pbs while
where appropriate. Otherwise, the executables in the algorithm's base directory
can be run with appropriate environment variables set (OMP\_NUM\_THREADS,
KMP\_AFFINITY, etc.).

### Inputs for executables
Running each program without arguments will print out the usage string.
Typically, each program takes as input serialized offset and edge list files
(IA, JA). Some of the algorithms accept additional arguments such as
delta-steps (SSSP), edge weights (VA files), and lists of source nodes.

IA, JA, and VA files can be created using the converter in matrix\_conversion/
and the initial .mtx file of the graph to be processed. The source lists can be
created from the .mtx file representing the graph sources to start from and the
sourcer.x executable, which is also compiled from matrix\_conversion/.
The .mtx files are available at /data/sparse.tamu.edu/gap/ on the Intel
devcloud, and the publically available, original versions are available 
[here](https://sparse.tamu.edu/GAP).