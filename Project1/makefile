CC = g++

CFLAGS = -std=c++11 -Wall

SOURCES = $(wildcard *.cpp) $(wildcard **/*.cpp)

OBJECTS = $(SOURCES:.c=.o)

TARGET = main

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

