FC = gfortran
FFLAGS = -O3

all: test

test: gamma_zeta_module.o test.o
	$(FC) $(FFLAGS) -o $@ $^

gamma_zeta_module.o: gamma_zeta_module.f90
	$(FC) $(FFLAGS) -c $<

test.o: test.f90 gamma_zeta_module.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod test

run: test
	./test

