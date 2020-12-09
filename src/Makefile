CC=g++
CFLAGS=-O3
LDFLAGS=-static

prog=a3m2msa fasta2pfam fastaNA fastaOneLine fastNf fixAlnX pfam2fasta RemoveNonQueryPosition trimBlastN


all: ${prog}


a3m2msa: a3m2msa.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastaNA: fastaNA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastaOneLine: fastaOneLine.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2pfam: fasta2pfam.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fastNf: fastNf.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fixAlnX: fixAlnX.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pfam2fasta: pfam2fasta.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

RemoveNonQueryPosition: RemoveNonQueryPosition.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

trimBlastN: trimBlastN.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${prog}
	cp ${prog} ../bin

clean:
	rm ${prog}
