CC = gcc
CFLAGS = -std=gnu11 -O -Wall 
LDLIBS = -lm

.PHONEY: default
default: out.txt sirPlot.png
	make sirPlot.png
	make harmonicPlot.png
	make threeBodyProblemPlot.png
	cat $<
	
sirPlot.png: sirPlot.gpi susceptible.txt infectious.txt removed.txt
	pyxplot $<

harmonicPlot.png: harmonicPlot.gpi harmonic.txt
	pyxplot $<

threeBodyProblemPlot.png: threeBodyProblem.gpi threeBodyProblem.txt
	pyxplot $<

out.txt: main
	./$< > $@

susceptible.txt: main
	./$< > $@

infectious.txt: main
	./$< > $@

removed.txt: main
	./$< > $@

harmonic.txt: main
	./$< > $@

main: main.o rungeKutta12.o
	$(CC) -o $@ $^ $(LDLIBS) 

main.o: main.c rungeKutta12.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

rungeKutta12.o: rungeKutta12.c rungeKutta12.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o *.txt *.png
