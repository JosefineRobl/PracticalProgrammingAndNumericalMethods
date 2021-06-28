CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt
	# Generates 'data.txt'
	make data.txt
	# Making the fit produces test to terminal, thus it is set first to make it easier to see out.txt text
	make timeDependence.png
	# Writes out the printed text from linearEquationSolver.c into the terminal
	cat out.txt
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

timeDependence.png: plot.gpi data.txt
	pyxplot $<

data.txt: linearEquationSolver Makefile
	# Passing argument to main in linearEquationSolver.c, thus creating data.txt. The argument is typed into the command line
	./$< "2"

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut* [Dd]ata* *.png
