CC = gcc
CFLAGS = -std=gnu11 -O -Wall 
LDLIBS = -lm
OBJFILES = main.o ode_rk12.o

default : sirPlot.png
	make sirPlot.png
	xdg-open sirPlot.png
	
sirPlot.png : plot.gpi susceptible.txt infectious.txt removed.txt
	pyxplot $<
	
susceptible.txt : main
	./$< > $@

infectious.txt: main
	./$< > $@

removed.txt: main
	./$< > $@

main : $(OBJFILES)
	$(CC) -o $@ $^ $(LDLIBS) 

main.o : main.c ode_rk12.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

ode_rk12.o : ode_rk12.c ode_rk12.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<



.PHONEY : clean 

clean : 
	@$(RM) main *.o *.txt *.png
