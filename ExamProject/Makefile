CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

EXE   	=  main5
OBJ  	=  main5.o integrator5.o # functions.o
HDRS    =  integrator5.h #functions.h
TXT     =  output.txt

.PHONEY: default
default: print

print: $(TXT)
	$(PIPE) $<

$(TXT) : $(EXE)
	./$< 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ)



#CC = gcc
#CFLAGS = -std=gnu11 -O -Wall
#LDLIBS = -lm

#.PHONEY: default
#default: out.txt
#	cat ./$<

#out.txt: main5
#	./$< > $@

#main5: main5.o integrator5.o
#	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

#main5.o: main5.c integrator5.h
#	$(CC) $(CFLAGS) -c $< -o $@

#integrator5.o: integrator5.c integrator5.h
#	$(CC) $(CFLAGS) -c $< -o $@

#.PHONEY: clean
#clean:
#	find -type f -executable -delete
#	$(RM) *.o [Oo]ut*
