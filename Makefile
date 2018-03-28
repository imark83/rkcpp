CC = g++
TARGETS = main.cpp fun.cpp rk.cpp buffer.cpp
OBJECTS = main.o fun.o rk.o buffer.o
LIBS =
CFLAGS = -O2 -g -Wall


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

test : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LIBS)


.PHONY: clean

clean:
	rm -f $(OBJECTS) test
