CFLAGS = -std=gnu11 -O -Wall
LDLIBS = lm

out.cmdline.txt: input.txt cmdline
	./read-cmdline $$(cat $<) > $@

input.txt: Makefile
	cat dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY: clean
clean:
	$(RM) input.txt out* cmdline

