# Program name and version
PROGRAM = StoCMet
VERSION = 1.0

# Directories
INST_DIR = $(HOME)/bin
SRC_DIR = .
BIN_DIR = ../bin
BUILD_DIR = ../build

# Compiler and flags
FC = gfortran


# Debug version    
#OPT = -g -fbounds-check -Wall -fcheck=all
#OPDE= d
# Optimized version
OPT = -O3
OPDE= o 

FCOPTS = -ffixed-line-length-132 -fcray-pointer -cpp -fdefault-real-8
FLAGS = -J$(BUILD_DIR) # Store .mod files in build directory

# Source files
SMODS = mod_Smolyak.f90

SFIX = mod_GaussHermite.f mod_GaussLegendre.f

SF90 = SCM.f90 lib_SCM.f90

# Object files
MODOBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SMODS))
FIXOBJ = $(patsubst %.f,$(BUILD_DIR)/%.o,$(SFIX))
F90OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SF90))

# Libraries
LDOPTO = 

# Default target
all: $(BUILD_DIR) program

# Rule to create the build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Linking the final executable
program: $(MODOBJ) $(FIXOBJ) $(F90OBJ) 
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -o $(BIN_DIR)/$(PROGRAM) $^ $(LDOPTO)

# Pattern rule to compile Fortran source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@

# Clean up generated files
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)/$(PROGRAM)


install:
	cp $(BIN_DIR)/$(PROGRAM) $(INST_DIR)/$(PROGRAM)-$(VERSION)-$(OPDE)	