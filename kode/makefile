# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags = c++ -std=c++0x -O3

MPIflags = mpic++ -std=c++11 -O3 
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG = main.exe
MPIPROG = MPImetropolis.exe
TEST_PROG = test_metropolis.exe

all : target tests MPItarget

target : ${PROG}
MPItarget : ${MPIPROG}
tests : ${TEST_PROG}

${PROG} : main.o metropolis.o
	${CPPflags} main.o metropolis.o ${LIB} -o ${PROG}

${TEST_PROG} : tests_main.o test_metropolis.o metropolis.o
	${CPPflags} tests_main.o test_metropolis.o metropolis.o ${LIB} -o ${TEST_PROG}

${MPIPROG} : MPImetropolis.o metropolis.o
	${MPIflags} MPImetropolis.o metropolis.o ${LIB} -o ${MPIPROG}

main.o : main.cpp
	${CPPflags} -c main.cpp

metropolis.o : metropolis.cpp
	${CPPflags} -c metropolis.cpp
	#${MPIflags} -c metropolis.cpp

MPImetropolis.o : MPImetropolis.cpp
	${MPIflags} -c MPImetropolis.cpp

tests_main.o : tests_main.cpp
	${CPPflags} -c tests_main.cpp

test_metropolis.o : test_metropolis.cpp
	${CPPflags} -c test_metropolis.cpp
