FC = gfortran
MYFFLAGS = -ffinite-math-only -funsafe-math-optimizations -ffast-math -funroll-loops  
FFLAGS = -Wall -Wextra -march=native $(MYFFLAGS)
MYLDFLAGS = -O3
LDFLAGS = $(MYLDFLAGS)

TARGETS = main
SOURCES = const.F90 integrators.F90 tidall.F90 run.F90 main.F90
OBJECTS = ${SOURCES:.f90=.o}

all: $(TARGETS)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

main: $(OBJECTS)
	$(FC) $(FFLAGS) -o tidal $^ $(LDFLAGS)

clean:
	@rm -rf *.o *.mod tidal

.PHONY: clean all

