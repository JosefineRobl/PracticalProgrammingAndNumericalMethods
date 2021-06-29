CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

.PHONEY: default
default: out.txt
	cat ./$<

out.txt: main4
	./$< > $@

main4: main4.o integrator4.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main4.o: main4.c integrator4.h
	$(CC) $(CFLAGS) -c $< -o $@

integrator4.o: integrator4.c integrator4.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
