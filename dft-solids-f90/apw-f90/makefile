



FC = gfortran

UTILS_DIR = $(HOME)/fortran-utils/src/


FCFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -I$(UTILS_DIR)
FLDFLAGS =  -llapack -L$(UTILS_DIR) -lfortran_utils


FSOURCES =   legendre_pols.f90 linalg.f90 quad.f90 utils.f90 apw.f90  double.f90  potential.f90


all: $(FSOURCES)
	$(FC)  $(FCFLAGS) $(FLDFLAGS) $? -o apw


 %.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<
