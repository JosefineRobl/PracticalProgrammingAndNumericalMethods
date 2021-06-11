CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt
	cat $<

out.txt: gslMatrix
	./$< > $@

gslMatrix: gslMatrix.o printMatrix.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

gslMatrix.o: gslMatrix.c printMatrix.h
	$(CC) $(CFSLGAS) -c $< -o $@

printMatrix.o: printMatrix.c printMatrix.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
