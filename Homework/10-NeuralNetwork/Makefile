CC = gcc
CFLAGS = -O -Wall -std=gnu11 $(shell gsl-config --cflags)
LDLIBS = -lm $(shell gsl-config --libs)

default: plot.png
	xdg-open $<

plot.png: plot.gpi dataFunctions.txt generatedPoints.txt
	pyxplot $<

out.txt: main
	./$< > $@

dataFunctions.txt: main
	./$<

generatedPoints.txt: main
	./$<

main: main.o artificialNeuralNetworks.o quasiNewton.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c artificialNeuralNetworks.h
	$(CC) $(CFLAGS) -c $< -o $@

artificialNeuralNetworks.o: artificialNeuralNetworks.c artificialNeuralNetworks.h quasiNewton.h
	$(CC) $(CFLAGS) -c $< -o $@

quasiNewton.o: quasiNewton.c quasiNewton.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY : clean 
clean :
	find -type f -executable -delete
	$(RM) *.txt *.o *.png
