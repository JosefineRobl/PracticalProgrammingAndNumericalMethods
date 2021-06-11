CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += -pthread
LDLIBS = -lm

.PHONEY: default
default: out.txt
	cat ./$<

out.txt: estimateOfPi
	./$< > $@

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) [Oo]ut*
