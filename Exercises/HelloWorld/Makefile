CC = gcc
CFLAGS = -O -std=gnu11 -march=native
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: helloWorld
	./helloWorld > out.txt

helloWorld: helloWorld.o
	$(CC) -o helloWorld helloWorld.o $(LDLIBS)

helloWorld.o: helloWorld.c
	$(CC) $(CFLAGS) -c helloWorld.c

clean: helloWorld.c
	$(RM) helloWorld.o helloWorld out.txt
