# Makefile: http://www.eng.hawaii.edu/Tutor/Make/

# Compiling flags
CFLAGS = -O3 -Wno-unused-result

# Linking flags
LFLAGS = -lm -lcpgplot -lpgplot -lX11 -fno-backslash -lpng

# Compilers
CC = gcc
F77 = gfortran

numint: numint.o vecmat.o force.o refsys.o rk.o jpleph.o
	$(CC) -o numint numint.o vecmat.o force.o refsys.o rk.o jpleph.o $(LFLAGS)

clean:
	rm -f *.o
	rm -f *~
