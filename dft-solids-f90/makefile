


FC = gfortran

FCFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace


FSOURCES =   apw.f90  double.f90 legendre_polynomial.f90


all: $(FSOURCES)
	$(FC)  $(FCFLAGS) $? -o dft


 %.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<
