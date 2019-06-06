CFLAGS += -fopenmp -O2

include $(BASILISK)/Makefile.defs

dirichlet.tst: CFLAGS += -DDIRICHLET=1
dirichlet.c: neumann_double.c
	ln -s neumann_double.c dirichlet.c
