FC = gfortran
MYFFLAGS = 
FFLAGS = -Wall -Wextra -march=native $(MYFFLAGS)
MYLDFLAGS = -O3
LDFLAGS = $(MYLDFLAGS)

TARGETS = main
SOURCES = auxs.f90 vars.f90 const.f90 funcs_f.f90 params_f.f90 derivates.f90 main.f90
OBJECTS = ${SOURCES:.f90=.o}

all: $(TARGETS)

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

main: $(OBJECTS)
	$(FC) $(FFLAGS) -o tidal $^ $(LDFLAGS)

clean:
	@rm -rf *.o *.mod tidal

.PHONY: clean all

