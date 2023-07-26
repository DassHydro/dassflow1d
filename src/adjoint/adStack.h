
#ifndef ADSTACK_LOADED
#define ADSTACK_LOADED 1

extern void pushNarray(char *x, unsigned int nbChars) ;

extern void popNarray(char *x, unsigned int nbChars) ;

extern void lookNarray(char *x, unsigned int nbChars) ;

extern void resetADLookStack() ;

extern void printbigbytes(long int nbblocks, long int blocksz, long int nbunits) ;

extern void printctraffic_() ;

extern void printtopplace_() ;

extern void printstackmax_() ;

extern void printlookingplace_() ;

extern void showrecentcstack_() ;

#endif
