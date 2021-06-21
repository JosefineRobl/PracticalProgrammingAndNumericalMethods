CC = gcc
CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm

default: out.txt
	cat ./$<

out.txt: main
	./$< > $@

main: main.o adaptiveIntegration2.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c adaptiveIntegration2.h
	$(CC) $(CFLAGS) -c $< -o $@

adaptiveIntegration2.o: adaptiveIntegration2.c adaptiveIntegration2.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
