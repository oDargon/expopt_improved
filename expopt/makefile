CC     = gcc
CFLAGS = -Wall -ggdb
SRC    = expopt.c inpctl.c energy.c opt.c amoeba.c mkscale.c qn.c \
         lnsrch.c gauss.c findmin.c guessian.c extend.c \
         chkhess.c jacobi.c powell.c chklindep.c log.c
OBJ    = ${SRC:.c=.o}
HDR    = expopt.h
BINDIR = /Users/widmark/bin

default: expopt2
all: expopt2 doc.pdf
expopt2: ${OBJ}
	${CC} -o expopt2 ${OBJ} -lm
doc.pdf: doc.tex fig.pdf
${OBJ}: ${HDR}
install: expopt2
	cp expopt2 ${BINDIR}/expopt2
clean:
	rm -fr expopt2 *.o *.lst *.tmp tmp.* tmp
	rm -fr *.aux *.dvi *.log *.pdf *.ps
.SUFFIXES:
.SUFFIXES: .eps .pdf .tex .c .o
.eps.pdf:
	epstopdf $*.eps
.tex.pdf:
	pdflatex $*
	pdflatex $*
.c.o:
	${CC} -c ${CFLAGS} $*.c
