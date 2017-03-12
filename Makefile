#
# Stampede - TACC
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = icpc
MPCC = mpicxx
OPENMP = -openmp
CFLAGS = -O3
LIBS =


TARGETS = serial pthreads openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
openmp: common.o grid.o openmp.o 
	$(CC) -o $@ $(OPENMP) openmp.o common.o grid.o $(LIBS)
mpi: mpicommon.o mpigrid.o mpi.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o mpicommon.o mpigrid.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h grid.h
	$(CC) -c $(OPENMP) $(CLFAGS)  openmp.cpp
serial.o: serial.cpp common.h 
	$(MPCC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp mpicommon.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp
	$(CC) -c $(CFLAGS) common.cpp
mpicommon.o:  mpicommon.cpp
	$(MPCC) -c $(CFLAGS) mpicommon.cpp
grid.o:	grid.cpp
	$(CC) -c $(CFLAGS) grid.cpp
mpigrid.o: mpigrid.cpp
	$(MPCC) -c $(CFLAGS) mpigrid.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
