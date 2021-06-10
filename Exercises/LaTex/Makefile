CC = gcc
CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

.PHONEY: default
default: report.pdf
	# Generates 'data.txt'
	make data.txt
	# Generate plot
	make plot.png
	# Generate pdf
	make $<
	# Open pdf
	xdg-open $<

report.pdf: report.tex plot.png
	pdflatex $<
	pdflatex $<

plot.png: plot.gpi data.txt
	pyxplot $<

data.txt: main
	./$< > $@

main: main.o exponential.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c exponential.h
	$(CC) $(CFLAGS) -c $< -o $@

exponential.o: exponential.c exponential.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean 
clean:
	find -type f -executable -delete
	$(RM) *.png *.pdf *.aux *.log *.o
