CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm
#CFLAGS += $(shell gsl-config --cflags)
#LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt error.png
	# Generates 'data.txt'
	make error.txt
	# Make the plot
	make error.png
	# Writes out the printed text from linearEquationSolver.c into the terminal
	cat out.txt
	# Opens the png image in the default image viewer on the computer
	xdg-open error.png

out.txt: monteCarloIntegration
	./$< > $@

monteCarloIntegration.o: monteCarloIntegration.c
	$(CC) $(CFLAGS) -c $< -o $@

monteCarloIntegration: monteCarloIntegration.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

error.png: plot.gpi error.txt
	pyxplot $<

error.txt: monteCarloIntegration Makefile
	# Passing argument to main in linearEquationSolver.c, thus creating data.txt. The argument is typed into the command line
	./$< "2"

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut* [Ee]rror* *.png
