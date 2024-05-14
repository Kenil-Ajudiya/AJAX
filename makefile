CUR_DIR=$(shell pwd)
#########################
OPENMP_INCL=
OPENMP_LIB=
ADDITIONAL_CFLAGS=
ADDITIONAL_FCFLAGS=
#########################

CC=g++
FC=gfortran

CFLAGS_DEBUG= -g -w -std=c++11

CFLAGS=-fPIC -c $(CFLAGS_DEBUG)  -I$(CUR_DIR)/include/  $(OPENMP_INCL) $(ADDITIONAL_CFLAGS) 

FCFLAGS=-L$(LD_LIBRARY_PATH) -fPIC -L$(CUR_DIR)/lib/ -lajax $(OPENMP_LIB) -fopenmp -L$(ADDITIONAL_FCFLAGS) -lgcc -lm -lc -lstdc++ -lX11 -lgfortran

all: libraries
	$(CC) $(CFLAGS) ajax.cpp -fopenmp
	$(FC) -o ajax ajax.o  $(FCFLAGS) 
	rm -rf *.o

libraries:
	make -C lib CC="$(CC)" CFLAGS=' $(CFLAGS)'

clean:
	make -C lib clean
	#make -C fileToSHM clean
	rm -rf ajax
