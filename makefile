# Makefile: http://www.eng.hawaii.edu/Tutor/Make/

# Compiling flags
CFLAGS = -O3 -Wno-unused-result

# Linking flags
LFLAGS = -lm -lcpgplot -lpgplot -lX11 -fno-backslash -lpng

# Compilers
CC = gcc
F77 = gfortran

tle2rv: tle2rv.o sgdp4.o satutl.o deep.o ferror.o refsys.o vecmat.o jpleph.o
	$(CC) -o tle2rv tle2rv.o sgdp4.o satutl.o deep.o ferror.o refsys.o jpleph.o vecmat.o -lm

numint: numint.o vecmat.o force.o refsys.o rk.o jpleph.o
	$(CC) -o numint numint.o vecmat.o force.o refsys.o rk.o jpleph.o $(LFLAGS)

clean:
	rm -f *.o
	rm -f *~
