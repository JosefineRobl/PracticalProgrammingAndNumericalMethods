CC = gcc
CFLAGS = -std=gnu11 -O -Wall 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt results.txt leastSquarePlot.png
	cat $<
	xdg-open leastSquarePlot.png

leastSquarePlot.png: results.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set title "Least square fit";\
	set xlabel "$$x$$";\
	set ylabel "$$y$$";\
	plot \
	"$<" index 0 using 1:2:3 with yerrorbars title "Data points",\
	"$<" index 1 using 1:2 with lines title "Fit",\
	"$<" index 2 using 1:2 with lines title "fit + $$c_0$$ and $$c_1$$ error",\
	"$<" index 3 using 1:2 with lines title "fit - $$c_0$$ and $$c_1$$ error"\
		' |pyxplot

results.txt: main
	./$< > $@

main: main.o leastSquareFunctions.o leastSquareFunctions.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c leastSquareFunctions.h
	$(CC) $(CFLAGS) -c $< -o $@

leastSquareFunctions.o: leastSquareFunctions.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY:
clean:
	find -type f -executable -delete
	$(RM) *.o *.txt *.png
