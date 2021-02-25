CXX = g++

CXXFLAGS = -std=c++11 -O1

SOURCES = $(wildcard *.cpp) $(wildcard **/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = main

$(TARGET) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^