CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt
	# Generates out.txt (this is due to the chosen way of placing the commands in this block, see reason below, since if 'make clean' used then 'data.txt' does not exist but shall be used in 'make timeDependence.png'
	make out.txt
	# Making the fit produces test to terminal, thus it is set first to make it easier to see out.txt text
	make timeDependence.png
	# Writes out the printed text from linearEquationSolver.c into the terminal
	cat out.txt
	#make timeDependence.png
	# Opens the png image in the default image viewer on the computer
	xdg-open timeDependence.png

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

timeDependence.png: plot.gpi data.txt Makefile
	pyxplot $<

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut* [Dd]ata* *.png
