CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

.PHONEY: default
default: out.txt
	cat ./$<

out.txt: main
	./$< > $@

main: main.o komplex.o # Using build-in linking rule

main.o komplex.o: komplex.h # Using build-in compiling rule

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
