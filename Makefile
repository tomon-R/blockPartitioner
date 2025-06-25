CXX = g++
CXXFLAGS = -std=c++17 -Wall

block_partitioner: block_partitioner.o
	$(CXX) $(CXXFLAGS) -o block_partitioner block_partitioner.o

block_partitioner.o: block_partitioner.cpp
	$(CXX) $(CXXFLAGS) -c block_partitioner.cpp

clean:
	rm -f block_partitioner block_partitioner.o