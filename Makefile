
CC = g++

CFLAGS = -g3 -O3 -ftree-vectorize -DHAVE_INLINE  -march=native -fopenmp -Igsl/include -std=c++11 

all: SDPR_admix score

SDPR_admix: parse_gen.o mcmc.o regress.o main.o
	${CC} ${CFLAGS} parse_gen.o mcmc.o regress.o main.o -Lgsl/lib/ -lgsl -lz -lgslcblas -lgomp -o SDPR_admix

main.o: main.cpp mcmc.h parse_gen.h mcmc.cpp parse_gen.cpp
	${CC} ${CFLAGS} -c main.cpp

regress.o: parse_gen.cpp parse_gen.h regress.cpp regress.h
	${CC} ${CFLAGS} -c regress.cpp

parse_gen.o: parse_gen.cpp parse_gen.h
	${CC} ${CFLAGS} -c parse_gen.cpp

mcmc.o: parse_gen.cpp parse_gen.h mcmc.cpp mcmc.h regress.cpp regress.h 
	${CC} ${CFLAGS} -c mcmc.cpp

score: score.o
	${CC} ${CFLAGS} score.o -lz -o score

score.o: score.cpp score.h
	${CC} ${CFLAGS} -c score.cpp

clean:
	rm -f *.o
