CC     = gcc
CFLAGS = -Wall -DDEBUG -DTRACE -g 
CFLAGS = -Wall
SRC    = onbas.c parse.c syntax.c rdbas.c prbas.c gram.c lowdin.c diag.c jacobi.c jacord.c
OBJ    = ${SRC:.c=.o}
LIB    = -lm
HDR    = onbas.h
BINDIR = /Users/widmark/bin

orth.exe: ${OBJ}
	${CC} -o orth.exe ${OBJ} ${LIB}
${OBJ}: ${HDR}
clean:
	rm -fr *.aux *.dvi *.log *.pdf *.ps *.bak *.o *.exe *.tmp tmp.* tmp
