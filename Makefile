FC = gfortran 
MYFFLAGS = -ffinite-math-only -funsafe-math-optimizations -ffast-math -funroll-loops 
FFLAGS = -Wall -Wextra -march=native -fcheck=all $(MYFFLAGS)
MYLDFLAGS = -O3
LDFLAGS = $(MYLDFLAGS)

TARGETS = main
FIXED_SOURCES = $(wildcard *.f) 
FREE_SOURCES = $(wildcard *.F90)
FIXED_OBJECTS = $(patsubst %.f,%.o,${FIXED_SOURCES})
FREE_OBJECTS = $(patsubst %.F90,%.o,${FREE_SOURCES})
DEP_FILE = my_project.dep

all: $(TARGETS) $(DEP_FILE) 

${FIXED_OBJECTS} : %.o : %.f
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

${FREE_OBJECTS} : %.o : %.F90
	$(FC) $(FFLAGS) -o $@ -c $< $(DEP_FILE) $(LDFLAGS)
	
%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

main: ${FREE_OBJECTS} ${FIXED_OBJECTS}
	$(FC) $(FFLAGS) -o tidal $^ $(LDFLAGS)

$(DEP_FILE):
	fortdepend -w -o $(DEP_FILE) -f *.F90

clean:
	@echo "rm -rf *.o *.mod tidal"
	@rm -rf *.o *.mod tidal

-include $(DEP_FILE)
depend: $(DEP_FILE)

.PHONY: clean

