CC = gcc
CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm

default: out.txt
	cat $<

out.txt: adaptiveIntegration
	./$< > $@

adaptiveIntegration: adaptiveIntegration.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)

adaptiveIntegration.o: adaptiveIntegration.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
