CC     = gcc
CFLAGS = -Wall -ggdb
SRC    = expopt.c inpctl.c energy.c opt.c amoeba.c mkscale.c qn.c qn2.c \
         lnsrch.c gauss.c sortxy.c findmin.c guessian.c extend.c \
         chkhess.c jacobi.c powell.c chklindep.c
OBJ    = ${SRC:.c=.o}
HDR    = expopt.h
BINDIR = /Users/widmark/bin

default: expopt
all: expopt doc.pdf
expopt: ${OBJ}
	${CC} -o expopt ${OBJ} -lm
doc.pdf: doc.tex fig.pdf
${OBJ}: ${HDR}
install: expopt
	cp expopt ${BINDIR}/expopt2
clean:
	rm -fr expopt *.o *.lst *.tmp tmp.* tmp
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
