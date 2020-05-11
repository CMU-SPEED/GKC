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
