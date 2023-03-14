CXX=gcc
CFLAG=-O3
LFLAG=-lpthread
EXE=pi.out
SRC=pi.c


OBJ=${SRC:.c=.o}

all: ${EXE}

${EXE}:${OBJ}
	${CXX} ${CFLAG} $^ -o $@ ${LFLAG}

%.o:%.c
	${CXX} ${CFLAG} -c $^ -o $@

.PHONY: clean

clean:
	rm -rf ${EXE} *.out ${OBJ}
