CFLAGS=-std=c++11 -O3 -I../common/ -Winline -Wall
PAR_FLAG=-fopenmp
ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG=-qopenmp
	CFLAGS+=-inline-forceinline -march=native -mavx2 -mavx512f
else # Assume g++
	CFLAGS+=-march=native -mavx2 -lm
endif

# Define Options:
# -DITERS=N to run N trials and return the average runtime
# -DDEBUG=N, for N in [0,3] for successively more verbose messaging
# -DBC_DUMP output betweenness centralities on stderr

all: bc bc_verify

bc: main.cpp bc.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CFLAGS} ${PAR_FLAG} $^ -o $@.exe

bc_verify: main.cpp bc.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CFLAGS} ${PAR_FLAG} -DVERIFY $^ -o $@.exe

clean: 
	rm -rf *.exe *.o
