#
#	Make file for program tcas.f
#

# location of fortran compiler

F77=gfortran

#---- start user specific make --------------------------------------------

tcas: tcas.o
	$(F77) -o tcas tcas.o

tcas.o: tcas.f
	$(F77) -c tcas.f
