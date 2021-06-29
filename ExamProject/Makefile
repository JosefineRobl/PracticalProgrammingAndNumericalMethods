CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

.PHONEY: default
default: out.txt
	cat ./$<

out.txt: main5
	./$< > $@

main5: main5.o integrator5.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main5.o: main5.c integrator5.h
	$(CC) $(CFLAGS) -c $< -o $@

integrator5.o: integrator5.c integrator5.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
