#Fortran compiler
FC = gfortran

#Compiler flags, -O3 optimizes the code to make it faster (-O0 is default, no optimization)
CFLAGS = -O3

#Source files
SOURCES = constants.f90 vectors.f90 simulation.f90

#Object files
OBJECTS = $(subst .f90,.o,$(SOURCES))

PROGRAM = simulation

$(PROGRAM): $(OBJECTS)
	$(FC) $(CFLAGS) -o $@ $^

%.o %.mod: %.f90
	$(FC) $(CFLAGS) -c -o $*.o $<
	@touch $@

vectors.o: constants.mod
Simulation.o: vectors.mod constants.mod


.PHONY: clean
clean: 
	rm -rf *.o *.mod $(PROGRAM)