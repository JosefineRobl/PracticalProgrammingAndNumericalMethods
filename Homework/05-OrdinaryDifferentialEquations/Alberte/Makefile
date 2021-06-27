CC = gcc
CFLAGS = std=gnu11 -O -Wall 
LDLIBS = -lm
OBJFILES = main.o ode_rk12.o
OUT = susceptible.txt infectious.txt removed.txt



default : plot.png
	
plot.png : plot.gpi $(OUT)
	pyxplot $<
	
$(OUT) : main
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
