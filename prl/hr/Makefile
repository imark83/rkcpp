CC = mpic++
# TARGETS = main.cpp fun.cpp rk.cpp work.cpp
OBJECTS = obj/main.o obj/fun.o obj/rk.o obj/work.o
CFLAGS = -O3 -Wall



aaa: clean $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)

obj/%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all: clean
	make aaa

lorenz: all
	cp aaa lorenz



.PHONY: clean

clean:
	rm -f $(OBJECTS) test
