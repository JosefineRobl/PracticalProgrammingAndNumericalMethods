CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: higgsMassPlot.png exerciseA.txt exerciseB.txt exerciseC.txt
	xdg-open $<

higgsMassPlot.png: higgsDataAllTogether.txt higgsMassData.txt Makefile
	echo '\
	set terminal png;\
	set key top right ;\
	set tics out;\
	set title "\textbf{Estimate of the mass of the Higgs boson}";\
	set xlabel "Energy E [GeV]";\
	set ylabel "Cross section $$\sigma(E)$$";\
	plot \
	"$<" using 1:2:3 with yerrorbars title "Data",\
	"higgsMassData.txt" index 0 using 1:2 with lines title "Fit with $$m = 125.97218$$";\
	' | pyxplot

exerciseA.txt: main
	./$< > $@

exerciseB.txt: main
	./$< > $@

exerciseC.txt: main
	./$< > $@

main.o: main.c quasiNewton.h simplex.h
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o quasiNewton.o simplex.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

quasiNewton.o: quasiNewton.c quasiNewton.h
	$(CC) $(CFLAGS) -c $< -o $@

simplex.o: simplex.c simplex.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) main *.o exercise* higgsMassData*
