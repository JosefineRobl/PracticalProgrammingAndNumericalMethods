CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: jacobiDiagonalisationTest.txt particleInABoxEigenvalues.txt particleInABox.png timePlot.png
	
particleInABox.png: particleInABox.gpi particleInABoxEigenfunctions.txt
	pyxplot $<

timePlot.png: timePlot.gpi diaganalisationTime.txt
	pyxplot $<

jacobiDiagonalisationTesting.txt: main
	./$< > $@

particleInABoxEigenvalues.txt: main
	./$< > $@

particleInABoxEigenfunctions.txt: main
	./$< > $@

diagonalisationTime.txt: main
	./$< > $@

main: main.o matrixDiagonalization.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c matrixDiagonalization.h
	$(CC) $(CFLAGS) -c $< -o $@

matrixDiagonalization.o: matrixDiagonalization.c matrixDiagonalization.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o *.txt *.png
