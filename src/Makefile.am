bin_PROGRAMS=phred
phred_SOURCES=phred.c \
	xmalloc.c \
	system.h \
	getopt.c \
	getopt1.c \
	getopt.h \
	lex.yy.c \
	jan-grammer.tab.c \
	Grid.cc \
	Material.cc \
	MaterialLib.cc 

EXTRA_DIST=phred.spec.in

yacc: jan-grammer.tab.c
	bison -d jan-grammer.y

lex: lex.yy.c
	flex jan.flex
