CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

TARGET = block_partitioner
SRC = block_partitioner.cpp
HDR = block_partitioner.hpp

all: $(TARGET)

$(TARGET): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
