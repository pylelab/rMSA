CC=g++
CFLAGS=-O3
LDFLAGS=-static

all: fastaNA catRNAcentral

fastaNA: fastaNA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

catRNAcentral: catRNAcentral.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
