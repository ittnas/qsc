Quantum solver collection - qsc
===============================
A collection of solvers designed for simulating quantum physics. Unfortunately, the documentation of the package is still lacking.

Installation
------------
First clone the repository ::

git clone https://github.com/ittnas/qsc/

After navigating to the repository build the library by running the makefile. For that you need a working version of lapack, blas and arpack. In Ubuntu, they can be simply installed as
'''
sudo apt install liblapack-dev
'''
Additionally, a fortran compiler is required. Install gfortan with ::
'''
sudo apt intall gfortran
'''
Additionally, sprng5 random number generator might be needed for the Monte-Carlo solvers. Its installation instructions can be found in /ext/sprng5.

In order to link to qsc you can use the following example makefile to build your simulation ::
'''
FC = gfortran
FCFLAGS = -Ofast -cpp -fopenmp
FCLIBS = -lqsc -llapack -lblas -lsprng -lstdc++
LIBPATH = -L$(QSCPATH)/bin -L$(QSCPATH)/ext/sprng5/lib
INCLPATH = -I$(QSCPATH)/include -I$(QSCPATH)/ext/sprng5/include

TARGET_PROGRAM = your_simulation

SRCS = $(TARGET_PROGRAM).f90
TARGET = $(TARGET_PROGRAM)

all: $(TARGET)

$(TARGET): $(SRCS)
	$(FC) $(FCFLAGS) $(SRCS) -o $@ $(INCLPATH) $(LIBPATH) $(FCLIBS)
'''
Here $(QSCPATH) is the environment variable containing the path to qsc root directory.