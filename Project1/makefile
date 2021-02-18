CC = g++

CFLAGS = -std=c++11 

SOURCES = $(wildcard *.cpp) $(wildcard **/*.cpp)

OBJECTS = $(SOURCES:.c=.o)

TARGET = main

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

