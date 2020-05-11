CXXFLAGS=-std=c++11 -O3 -I../common -Winline -march=native -mavx2
PAR_FLAG=-fopenmp
ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG=-qopenmp
	CXXFLAGS+=-inline-forceinline -mavx512f
endif

all: pagerank

pagerank: pagerank.c pagerank.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} $^ -o $@.exe

clean:
	rm -f *.o *.exe
