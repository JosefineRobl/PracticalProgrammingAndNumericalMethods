CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: exerciseAPlot.png exerciseBPlot.png exerciseCPlot.png

exerciseAPlot.png: exerciseAPlot.gpi xyPoints.txt exerciseA.txt 
	pyxplot $<

qplot.png: exerciseBPlot.gpi xyPoints.txt exerciseB.txt
	pyxplot $<

cplot.png: exerciseCPlot.gpi xyPoints.txt exerciseC.txt
	pyxplot $<

exerciseA.txt: main
	./$< > $@

exerciseB.txt: main
	./$< > $@

exerciseC.txt: main
	./$< > $@

main.o: main.c linearSpline.h cubicSpline.h quadraticSpline.h
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o linearSpline.o cubicSpline.o quadraticSpline.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

linearSpline.o: linearSpline.c linearSpline.h
	$(CC) $(CFLAGS) -c $< -o $@

cubicSpline.o: cubicSpline.c cubicSpline.h
	$(CC) $(CFLAGS) -c $< -o $@

quadraticSpline.o: quadraticSpline.c quadraticSpline.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	find -type f -executable -delete
	$(RM) *.o *.txt
