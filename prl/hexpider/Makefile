CC = mpic++
# TARGETS = main.cpp fun.cpp rk.cpp work.cpp
OBJECTS = obj/main.o obj/fun.o obj/rk.o
CFLAGS = -O3 -Wall -std=c++11



aaa: clean $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)

obj/%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all: clean
	make aaa



.PHONY: clean

clean:
	rm -f $(OBJECTS) test
