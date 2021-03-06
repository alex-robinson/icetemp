.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = .
objdir = .obj
bindir = .

# Command-line options at make call
debug ?= 0

## LIB CONFIGURATION ## 
INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

FC = gfortran

FFLAGS  = -ffree-line-length-none -I$(objdir) -J$(objdir) $(INC_NC)
LFLAGS  = $(LIB_NC)
DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg


# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of rules and source files
##
###############################################

$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/defs.o: $(srcdir)/defs.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_tridiagonal.o: $(srcdir)/solver_tridiagonal.f90 $(objdir)/defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/thermodynamics.o : $(srcdir)/thermodynamics.f90 $(objdir)/defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_column_interp.o : $(srcdir)/ice_column_interp.f90 $(objdir)/defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_temperature.o : $(srcdir)/ice_temperature.f90 $(objdir)/defs.o \
					  $(objdir)/solver_tridiagonal.o $(objdir)/thermodynamics.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_enthalpy.o : $(srcdir)/ice_enthalpy.f90 $(objdir)/defs.o \
					  $(objdir)/solver_tridiagonal.o $(objdir)/thermodynamics.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ice_age.o : $(srcdir)/ice_age.f90 $(objdir)/defs.o \
					  $(objdir)/solver_tridiagonal.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
##
## Compilation of complete programs
##
###############################################

test_icetemp : $(objdir)/ncio.o $(objdir)/defs.o $(objdir)/solver_tridiagonal.o \
				$(objdir)/thermodynamics.o $(objdir)/ice_temperature.o $(objdir)/ice_enthalpy.o \
				$(objdir)/ice_age.o $(objdir)/ice_column_interp.o
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/test_icetemp.x test_icetemp.f90 \
			$(LFLAGS) $^
		@echo " "
		@echo "    test_icetemp.x is ready."
		@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make test_icetemp : compiles test_icetemp.x, stand-alone icetemp test program."
	@echo " make clean     : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
