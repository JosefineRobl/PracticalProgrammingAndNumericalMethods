CC = gcc
CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm

default: out.txt
	cat ./$<

out.txt: main
	./$< > $@

main: main.o integratorTridivision.o integratorBidivision.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c integratorTridivision.h integratorBidivision.h
	$(CC) $(CFLAGS) -c $< -o $@

integratorTridivision.o: integratorTridivision.c integratorTridivision.h
	$(CC) $(CFLAGS) -c $< -o $@

integratorBidivision.o: integratorBidivision.c integratorBidivision.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut*
