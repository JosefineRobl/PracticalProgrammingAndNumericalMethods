CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: boundStatesHydrogen.png exerciseA.txt exerciseB.txt
	xdg-open $<	

exerciseA.txt: main
	./$< > $@

exerciseB.txt: main
	./$<

boundStatesHydrogen.txt: main
	./$<

boundStatesHydrogen.png: boundStatesHydrogen.gpi boundStatesHydrogen.txt
	pyxplot $<

main: main.o qrGramSchmidt.o rungeKutta23.o rootFinding.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

linearEquationSolver.o: linearEquationSolver.c qrGramSchmidt.h rungeKutta23.h rootFinding.h
	$(CC) $(CFLAGS) -c $< -o $@

qrGramSchmidt.o: qrGramSchmidt.c qrGramSchmidt.h
	$(CC) $(CFLAGS) -c $< -o $@

rungeKutta23.o: rungeKutta23.c rungeKutta23.h
	$(CC) $(CFLAGS) -c $< -o $@

rootFinding.o: rootFinding.c rootFinding.h qrGramSchmidt.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -deleta
	$(RM) *.o *.txt *.png
