all: converter sourcer

sourcer: mmio.c read_sources.cpp
	g++ -g --std=c++11 -O2 $^ -o $@.x

converter: mmio.c read_mmio.c
	g++ -g --std=c++11 -O2 -DDEBUG $^ -o $@.x

clean:
	rm -rf *.x *.o
