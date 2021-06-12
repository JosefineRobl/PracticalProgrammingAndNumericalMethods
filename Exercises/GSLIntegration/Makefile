CC = gcc
CFLAGS = -std=gnu11 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.txt errfun.png #bessel.png
	cat ./$<
	make dataB.txt
	#make dataC.txt
	make errfun.png
	#make bessel.png
	xdg-open errfun.png
	#xdg-open bessel.png

out.txt: gslIntegration
	./$< > $@

dataB.txt: gslIntegration

#dataC.txt: gslIntegration

gslIntegration: gslIntegration.o
	$(CC) $(LDFLAGS) $< -o $@ $(LDLIBS)

gslIntegration.o: gslIntegration.c
	$(CC) $(CFLAGS) -c $< -o $@

errfun.png: errfun.gpi dataB.txt
	pyxplot $<

#bessel.png: bessel.gpi dataC.txt
#	pyxplot $<

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o [Oo]ut* [Dd]ata* *.png
