CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

.PHONEY: default
default: out.txt
	cat ./$<

out.txt: main
	./$< > $@

main: main.o integrator.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c integrator.h
	$(CC) $(CFLAGS) -c $< -o $@

integrator.o: integrator.c integrator.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
