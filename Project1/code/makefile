CXX = g++

CXXFLAGS = -std=c++11 -O3 -fopenmp

SOURCES = $(wildcard *.cpp) $(wildcard **/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = main

$(TARGET) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^