CFLAGS=-std=c++11 -O3 -I../common/ -Winline -Wall
PAR_FLAG=-fopenmp
ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG=-qopenmp
	CFLAGS+=-inline-forceinline -march=native -mavx2 -mavx512f
else # Assume g++
	CFLAGS+=-march=native -mavx2 -lm
endif


# Define Options:
# -DCHECK_DEPTHS to verify parents and depths against seqential
# -DITERS=N to run N trials and return the average runtime

all: bfs bfs_verify

bfs: bfs.cpp bfs.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CFLAGS} ${PAR_FLAG} -DITERS=1 $^ -o $@.exe

bfs_verify: bfs.cpp bfs.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CFLAGS} ${PAR_FLAG} -DITERS=1 -DCHECK_DEPTHS $^ -o $@.exe

clean: 
	rm -rf *.exe *.o
