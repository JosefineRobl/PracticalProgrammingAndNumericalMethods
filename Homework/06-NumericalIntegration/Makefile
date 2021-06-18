CC = gcc
CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm

default: out.txt
	cat ./$<

out.txt: main
	./$< > $@

main: main.o adaptiveIntegration.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c adaptiveIntegration.h
	$(CC) $(CFLAGS) -c $< -o $@

adaptiveIntegration.o: adaptiveIntegration.c adaptiveIntegration.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
