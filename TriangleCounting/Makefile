CXXFLAGS=-std=c++11 -O3 -I../common -Winline
PAR_FLAG=-fopenmp
ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG=-qopenmp
else # Assume g++
	PAR_FLAG=-fopenmp
	CXXFLAGS+=-march=native
endif

# Additional options:
# -DITERS=1

all: tc tc_verify

tc: main.c tc.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} $^ -o $@.exe

tc_verify: main.c tc.a ../common/graph.cpp ../common/utils.cpp
	${CXX} ${CXXFLAGS} ${PAR_FLAG} -DVERIFY $^ -o $@.exe

clean: 
	rm -rf *.o *.exe
