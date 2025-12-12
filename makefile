all: aout

aout: mod_spectra_cmrt.o spectra_cmrt.o
	ifort -mkl mod_spectra_cmrt.f90 spectra_cmrt.f90 -o aout

%.o: %.f90
	ifort -c $<

clean:
	rm *.o aout

