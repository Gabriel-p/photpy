Para compilar:


g77 -c (soubroutine).f (subroutine).o :

g77 -c daomaster.f daomaster.o

g77 -c mathsubs.f mathsubs.o

g77 -c iosubs.f iosubs.o

g77 -c unxsubs.f unxsubs.o




g77 -o (ejecutable) (programa).f (subroutine1).o (subroutine2).o ..... :

g77 -o daomaster daomaster.o iosubs.o mathsubs.o unxsubs.o

