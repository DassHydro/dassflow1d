
SCOTCHDIR  = $(libdir)/../../scotch_5.1.12_esmumps
ISCOTCH    = -I$(SCOTCHDIR)/include

LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

ORDERINGSF  = -Dpord -Dscotch
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LPORD) $(LSCOTCH)
IORDERINGSF = $(IPORD) $(ISCOTCH)
IORDERINGSC = $(IPORD) $(ISCOTCH)

PLAT    =
LIBEXT  = .a
OUTC    = -o 
OUTF    = -o 
RM = /bin/rm -f
CC = gcc
FC = gfortran
FL = gfortran
AR = ar vr 

RANLIB  = echo
INCSEQ = -I$(topdir)/libseq
LIBSEQ  =  -L$(topdir)/libseq -lmpiseq
LIBBLAS = -lblas
LIBOTHERS = -lpthread
CDEFS   = -DAdd_

OPTF    = -O3  -DALLOW_NON_INIT -fPIC -fallow-argument-mismatch
OPTL    = -O3 -fPIC
OPTC    = -O3 -fPIC

INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
