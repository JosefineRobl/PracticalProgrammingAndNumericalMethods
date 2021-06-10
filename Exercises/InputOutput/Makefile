CC = gcc
CFLAGS = -std=gnu11 -O -Wall
LDLIBS = -lm

.PHONEY: default
default: out.txt
	cat ./$<

out.txt: out.cmdline.txt out.stdio.txt out.file.txt
	cat /dev/null > $@
	cat $^ >> $@
#	cat $< > $@
#	cat out.stdio.txt >> $@
#	cat out.file.txt >> $@

out.cmdline.txt: cmdline Makefile
	./cmdline 1 2 3 4 5 6 7 8 9 10 > $@

cmdline: cmdline.o printExercise.o formatPrint.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

cmdline.o: cmdline.c printExercise.h
	$(CC) $(CFLAGS) -c $< -o $@

out.stdio.txt: stdio Makefile
	@echo 0 10 20 30 40 50 60 70 80 90 | ./stdio > $@

stdio: stdio.o printExercise.o formatPrint.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

stdio.o: stdio.c printExercise.h
	$(CC) $(CFLAGS) -c $< -o $@

out.file.txt: input.txt file
	./file "$<" > $@

file: file.o printExercise.o formatPrint.o
	$(CC) $(LSFLAGS) $^ -o $@ $(LDLIBS)

file.o: file.c printExercise.h
	$(CC) $(CFLAGS) -c $< -o $@

printExercise.o: printExercise.c printExercise.h
	$(CC) $(CFLAGS) -c $< -o $@

formatPrint.o: formatPrint.c formatPrint.h
	$(CC) $(CFLAGS) -c $< -o $@

input.txt: Makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY: clean
clean:
	find -type f -executable -delete
	$(RM) *.o input.txt [Oo]ut*

