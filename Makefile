# Fortran Compiler
FC = gfortran

# Library Paths and Linking Flags
LDFLAGS = -L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/lapack/lib -lopenblas -llapack

# Source Files
SOURCES = Heis.f90 Det.f90 Factorial.f90 ReadInput.f90 Integ_to_Bit.f90 ValCoupl.f90 Diasym.f90 Identify_Ms.f90 QEqual.f90

# Object Files
OBJECTS = $(SOURCES:.f90=.o)

# Phony Targets
.PHONY: all clean

# Default Target
all: Heis

# Linking the Executable
Heis: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(LDFLAGS)

# Compiling Source Files to Object Files
%.o: %.f90
	$(FC) -c $<

# Clean Target
clean:
	rm -rf *.o Heis

