CFLAGS +=  -O2

include $(BASILISK)/Makefile.defs

dirichlet.tst: CFLAGS += -DDIRICHLET=1 -DMTRACE=3
dirichlet.c: neumann_double.c
	ln -s neumann_double.c dirichlet.c
