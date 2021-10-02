FC = gfortran
MYFFLAGS = 
FFLAGS = -Wall -Wextra -march=native $(MYFFLAGS)
MYLDFLAGS = -O3
LDFLAGS = $(MYLDFLAGS)

TARGETS = main
SOURCES = const.f90 f_funcs.f90 tidal_vars.f90 tidal_derivates.f90 integrators.f90 run.f90 main.f90
OBJECTS = ${SOURCES:.f90=.o}

all: $(TARGETS)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

main: $(OBJECTS)
	$(FC) $(FFLAGS) -o tidal $^ $(LDFLAGS)

clean:
	@rm -rf *.o *.mod tidal

.PHONY: clean all

