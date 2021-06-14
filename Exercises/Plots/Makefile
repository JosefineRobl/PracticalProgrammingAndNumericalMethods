CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: erf.png gamma.png
	make dataErf.txt
	make dataGamma.txt
	make erf.png
	make gamma.png

erf.png: erf.gpi dataErf.txt tabulatedDataErf.txt
	pyxplot $<

gamma.png: gamma.gpi dataGamma.txt tabulatedDataGamma.txt
	pyxplot $<

dataErf.txt: main
	./$< "1" > $@

dataGamma.txt: main
	./$< "2" > $@

main: main.o erfApproximation.o gammaApproximation.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c erfApproximation.h gammaApproximation.h
	$(CC) $(CFLAGS) -c $< -o $@

erfApproximation.o: erfApproximation.c erfApproximation.h
	$(CC) $(CFLAGS) -c $< -o $@

gammaApproximation.o: gammaApproximation.c gammaApproximation.h
	$(CC) $(CFLAGS) -c $< -o $@

dataErfTabulated.txt:

dataGammaTabulated.txt:

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o *.png [Dd]ata*
