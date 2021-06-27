CC = gcc
CFLAGS = -std=gnu11 -O -Wall 
LDLIBS = -lm
OUT = susceptible.txt infectious.txt removed.txt out.txt

.PHONEY: default
default: out.txt
	make $(OUT)
	make sirPlot.png
	cat $<

sirPlot.png: sirPlot.gpi $(OUT)
	pyxplot $<

$(OUT): main
	./$< > $@

main: main.o rungeKutta12.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c rungeKutta12.h
	$(CC) $(CFLAGS) -c $< -o $@

rungeKutta12: rungeKutta12.c rungeKutta12.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.txt *.png
