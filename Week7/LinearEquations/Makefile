CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: linearEquationSolver
	./$< > $@

linearEquationSolver.o: linearEquationSolver.c
	$(CC) $(CFLAGS) -c $< -o $@

linearEquationSolver: linearEquationSolver.o
	$(CC) $(LDFLAGS) linearEquationSolver.o -o linearEquationSolver $(LDLIBS)

clean:
	$(RM) *.o out* linearEquationSolver
