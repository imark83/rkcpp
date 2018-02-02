CC = g++
TARGETS = main.cpp fun.cpp ../rk.cpp
OBJECTS = main.o fun.o ../rk.o
LIBS = -fopenmp
CFLAGS = -O2 -g -Wall -fopenmp


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

aaa : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LIBS)

fli:
	make fli



.PHONY: clean

clean:
	rm -f $(OBJECTS)
