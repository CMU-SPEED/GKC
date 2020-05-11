CXXFLAGS=-std=c++11 -O3 -march=native -mavx2 -I../common/ -Winline
PAR_FLAG=-fopenmp
ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG=-qopenmp
	CXXFLAGS+=-inline-forceinline -mavx512f 
else # Assume g++
	PAR_FLAG=-fopenmp
endif

# Additional options:
# -DITERS=1

all: conn_comps conn_comps_verify

conn_comps: main.cpp conn_comps.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} $^ -o $@.exe

conn_comps_verify: main.cpp conn_comps.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} -DVERIFY $^ -o $@.exe

clean: 
	rm -rf *.o *.exe
