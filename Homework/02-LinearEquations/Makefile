CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: linearEquationSolver
	./$< > $@

linearEquationSolver.o: linearEquationSolver.c qrGramSchmidt.h printMatrix.h
	$(CC) $(CFLAGS) -c $< -o $@

qrGramSchmidt.o: qrGramSchmidt.c qrGramSchmidt.h
	$(CC) $(CFLAGS) -c $< -o $@

printMatrix.o: printMatrix.c printMatrix.h
	$(CC) $(CFLAGS) -c $< -o $@

linearEquationSolver: linearEquationSolver.o qrGramSchmidt.o printMatrix.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
