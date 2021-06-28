CC = gcc
CFLAGS = -O -Wall -std=gnu11 $(shell gsl-config --cflags)
LDLIBS = -lm $(shell gsl-config --libs)

default: plot.png
	xdg-open $<

plot.png: plot.gpi dataFunctions.txt points.txt
	pyxplot $<

dataFunctions.txt generatedPoints.txt log.txt: main
	./$< > data.txt points.txt 2> out.txt

main: main.o artificialNeuralNetwork.o quasiNewton.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c artificialNeauralNetwork.h
	$(CC) $(CFLAGS) -c $< -o $@

artificialNeauralNetwork.o: artificialNeauralNetwork.c artificialNeauralNetwork.h quasiNewton.h
	$(CC) $(CFLAGS) -c $< -o $@

quasiNewton.o: quasiNewton.c quasiNewton.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY : clean 
clean :
	find -type f -executable -delete
	$(RM) *.txt *.o *.png
