# Graph Kernel Collection
#
# Copyright 2020 Carnegie Mellon University.
#
# NO WARRANTY. THIS CARNEGIE MELLON UNIVERSITY AND SOFTWARE ENGINEERING
# INSTITUTE MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. CARNEGIE MELLON 
# UNIVERSITY MAKES NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, 
# AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR 
# PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF 
# THE MATERIAL. CARNEGIE MELLON UNIVERSITY DOES NOT MAKE ANY WARRANTY OF ANY
# KIND WITH RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT 
# INFRINGEMENT.
#
# Released under a BSD (SEI)-style license, please see license.txt or
# contact permission@sei.cmu.edu for full terms.
#
# [DISTRIBUTION STATEMENT A] This material has been approved for public
# release and unlimited distribution.  Please see Copyright notice for 
# non-US Government use and distribution.
#
# This Software includes and/or makes use of the following Third-Party
# Software subject to its own license:
#
# 1. Matrix Market Loader code (https://math.nist.gov/MatrixMarket/mmio-c.html).
#
#      The code made publicly available at nist.gov is not marked with a 
#      copyright notice and is therefore believed pursuant to section 105 of 
#      the Copyright Act, to not be entitled to domestic copyright protection 
#      under U.S. law and is therefore in the public domain.  Accordingly, it 
#      is believed that no license is required for its use.
#
# This Software may include certain portions of copyrighted code that is 
# initially being released only in binary form for validation and evaluation
# purposes. It is expected that source code will be released as open source at
# a future date. 
#
# DM20-0375

CXXFLAGS=-std=c++11 -O3 -I../common -Winline -march=native -mavx2
PAR_FLAG=-fopenmp
ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG=-qopenmp
	CXXFLAGS+=-inline-forceinline -mavx512f
endif

# Additional options:
# -DITERS=1

DATAPATH=/sharedstorage/markb1/GAP_data/processed
all: sssp sssp_verify

# Option to dump verified values: -DDUMP_DISTS
# These will be dumped to stderr
# Option to control number of sources processed:
# -DITERS=<INT>
# -DVALIDATE enables checking node distances against a sequential dijkstra
#  implementation.

sssp: sssp.cpp sssp.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} $^ -o $@.exe

sssp_verify: sssp.cpp sssp.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} -DVALIDATE $^ -o $@.exe

clean: 
	rm -rf *.o *.exe
