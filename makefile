FC = gfortran
#FCFLAGS = -c -Ofast -fopenmp -cpp -fbounds-check -Wall -finit-real=nan
#FCFLAGS = -c -Ofast -fopenmp -cpp -fbounds-check -finit-real=nan
FCFLAGS = -c -Ofast -fopenmp -cpp
#FCFLAGS = -c -O -fopenmp -cpp -g
FCLIBS = -llapack -lblas -larpack
EXTLIBS = -L
EXTINCL = -I
#EXTLIBS = -Lext/SOFTWARE -lSparseBLAS_GNU
#EXTINCL = -Iext/SOFTWARE
#SPRNGLIBS = -Lext/sprng5/SRC -lsprng -lstdc++
SPRNGLIBS = -Lext/sprng5/lib -lsprng -lstdc++
SPRNGINCL = -Iext/sprng5/include

LD = ar
LDFLAGS = -crs

SRCDIR = src
MODDIR = include
OBJDIR = obj
BINDIR = bin

SRCS = utils.f90 constants.f90 random.f90 linfuncs.f90 operators.f90 envelopes.f90 signals.f90 hamiltonians.f90 sclass.f90 fclass.f90 equation_solvers.f90 ode_solvers.f90 sde_solvers.f90 ode_s_solvers.f90

SOURCES  := $(wildcard $(addprefix $(SRCDIR)/,$(SRCS)))
#INCLUDES := $(wildcard $(SRCDIR)/*.h)
#OBJECTS  := $(SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)

#OBJS = $(SRCS:.f90=.o)
LIBNAME = libqsc

$(BINDIR)/$(LIBNAME): $(OBJECTS)
	@$(LD) $(LDFLAGS) $@.a $(OBJECTS)
	@echo "Linking complete!"

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@$(FC) $(FCFLAGS) $< $(FCLIBS) -o $@ -J$(MODDIR)/ -I$(MODDIR)/ $(EXTLIBS) $(EXTINCL) $(SPRNGLIBS) $(SPRNGINCL)
	@echo "Compiled "$<" successfully!"



#all: qsl
#all: $(SRCS) $(LIBNAME)
#$(LIBNAME): $(OBJS) 
#	$(LD) $(LDFLAGS) $@.a $(OBJS)
#%.o: %.f90
#	$(FC) $(FCFLAGS) $< $(FCLIBS)
