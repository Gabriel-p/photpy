F77=gfortran
FFLAGS = -c -O3
LFLAGS = -O3

# RULES:
.SUFFIXES: .o .f
.f.o:
	$(F77) $(FFLAGS) $<

daomaster: daomaster.o iosubs.o mathsubs.o unxsubs.o
	f77 $(LFLAGS) -o daomaster daomaster.o iosubs.o mathsubs.o unxsubs.o

all:
	make daomaster

clean:
	rm *.o

