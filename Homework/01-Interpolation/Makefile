CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: exerciseAPlot.png exerciseBPlot.png exerciseCPlot.png

exerciseAPlot.png: exerciseAPlot.gpi xyPoints.txt exerciseA.txt 
	pyxplot $<

exerciseBPlot.png: exerciseBPlot.gpi xyPoints.txt exerciseB.txt
	pyxplot $<

exerciseCPlot.png: exerciseCPlot.gpi xyPoints.txt exerciseC.txt
	pyxplot $<

exerciseA.txt: main
	./$< > $@

exerciseB.txt: main
	./$< > $@

exerciseC.txt: main
	./$< > $@

xyPoints.txt: main
	./$< > $@

main.o: main.c linearSpline.h cubicSpline.h quadraticSpline.h binarySearch.h
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o linearSpline.o cubicSpline.o quadraticSpline.o binarySearch.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

linearSpline.o: linearSpline.c linearSpline.h binarySearch.h
	$(CC) $(CFLAGS) -c $< -o $@

cubicSpline.o: cubicSpline.c cubicSpline.h binarySearch.h
	$(CC) $(CFLAGS) -c $< -o $@

quadraticSpline.o: quadraticSpline.c quadraticSpline.h binarySearch.h
	$(CC) $(CFLAGS) -c $< -o $@

binarySearch.o: binarySearch.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o *.txt
