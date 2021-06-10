FLAGS = -o -w -std=gnu99
CFLAGS += $(shell gsl-config --cflags)
#CFLAGS += 'gsl-config --cflags' # One needs to use the shell command since the ' ' commands does not work anymore.
LDLIBS += $(shell gsl-config --libs)
#LDLIBS += $(gsl-config --libs) # One needs to add shell to use the terminal

all: output ; cat $<
#all: output ; cat $@ # $@ is the target while $< is the first prerequisite
output: main ; ./$< > $@
#	output: main : $@ > $< # Shall not be indented since it is the "header" for the rule and not whats in it. Furthermore shall : -> ; since semi-colon is used instead of indentation to show when the rule header ends and the body begins. The correct way is: target:prerequisite;rule. Furthermore the rule shall have switched $< and $@, since $@ is the target to which the prerequisite ($<) shall be "pushed" into (>). 
clean: ; $(RM) main output
#clean: ; $RM main output # There shall be a paranthesis around RM, since RM together is the variable.
