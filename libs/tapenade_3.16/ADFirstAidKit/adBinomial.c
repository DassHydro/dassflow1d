/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2021 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 */

#include <stdio.h>
#include <stdlib.h> //exit
#include "adBinomial.h"

// Max number of nested binomial checkpoint loops:
#define _ADBINOMIALMAXNEST 5
// Max number of snapshots in a binomial checkpoint loop:
#define _ADBINOMIALMAXSNP 99

static int stack1[3*_ADBINOMIALMAXNEST] ;
static int stack2[3*_ADBINOMIALMAXSNP] ;
static int is1 = -1 ;
static int is2 = -1 ;

#define BOTTOM stack1[3*is1]
#define MAXSNP stack1[3*is1+1]
#define OFFSET stack1[3*is1+2]
#define LENGTH stack2[3*is2]
#define CURPOS stack2[3*is2+1]
#define CKPCUT stack2[3*is2+2]

char *actionText[] = {"","PUSHSNAP","LOOKSNAP","POPSNAP","ADVANCE","FIRSTTURN","TURN"} ;

void adBinomial_init(int length, int nbSnap, int firstStep) {
  if (length<=0) {
    printf("Error: Cannot reverse a sequence of length %i\n",length) ;
    exit(0) ;
  } else if (nbSnap==0 && length>=2) {
    printf("Error: Cannot reverse a sequence of length %i with no snapshot\n",length) ;
    exit(0) ;
  } else if (is1>=_ADBINOMIALMAXNEST || is2+nbSnap+1>=_ADBINOMIALMAXSNP) {
    printf("Error: Binomial-Checkpointing memory exceeded !\n") ;
    exit(0) ;
  } else {
    is1 = is1+1 ;
    is2 = is2+1 ;
    BOTTOM = is2 ;
    MAXSNP = nbSnap ;
    OFFSET = firstStep ;
    LENGTH= length ;
    CURPOS = 0 ;
    CKPCUT = 0 ;
  }
}

void adBinomial_setCut() {
  int length = LENGTH ;
  int nbSnap = MAXSNP-(is2-BOTTOM)+1 ;
  if (length<=1) {
    CKPCUT = 0 ;
  } else if (nbSnap==1) {
    CKPCUT = length-1 ;
  } else {
    int minRecomp = 1 ;
    int eta = nbSnap+1 ;
    while (eta<length) {
      minRecomp++ ;
      eta = (eta*(nbSnap+minRecomp))/minRecomp ;
    }
    CKPCUT = (length*minRecomp)/(minRecomp+nbSnap) ;
    if (CKPCUT==0)
      CKPCUT=1 ;
    else if (CKPCUT>=length)
      CKPCUT=length-1 ;
  }
}

int adBinomial_next(int *action, int *step) {
  int i ;

//  Begin "only for debug" part:
//  printf("\n") ;
//  printf("STACK1:\n") ;
//  for (i=is1 ; i>=0 ; i--)
//    printf("%i snapshots, stack2 bottom:%i (offset:%i)\n",
//           stack1[3*i+1],stack1[3*i],stack1[3*i+2]) ;
//  printf("-------------------\n") ;
//  printf("STACK2:\n") ;
//  for (i=is2 ; i>=0 ; i--)
//    printf("%i: R( ,%i) %i/%i\n",
//           i,stack2[3*i],stack2[3*i+1],stack2[3*i+2]) ;
//  //    printf("%i: at %i, nextckp %i, end %i\n",
//  //           i,stack2[3*i+1],stack2[3*i+2],stack2[3*i]) ;
//  printf("-------------------\n") ;
//  End "only for debug" part.

  if (LENGTH<=0 && is2==BOTTOM) {
    *step = -1 ;
    *action = -1 ;
    is1-- ;
    is2-- ;
// printf("NEXT BINOMIAL ACTION: -1\n") ;
    return 0 ;
  } else {
    *step = 1 ;
    for (i=BOTTOM+1 ; i<=is2 ; i++) *step += stack2[3*i+1] ;
    if (CURPOS==-1) {
      *action = (LENGTH==1)?POPSNAP:LOOKSNAP ;
      CURPOS = 0 ;
    } else {
      if (CURPOS==LENGTH-1) {
        *action = (*step==stack2[3*BOTTOM])?FIRSTTURN:TURN ;
        if (CURPOS==0 && is2>BOTTOM) {
          is2-- ;
          LENGTH = CKPCUT ;
        } else {
          LENGTH = CURPOS ;
        }
        CURPOS = -1 ;
        adBinomial_setCut() ;
      } else if (CURPOS==CKPCUT) {
        int remainingLength = LENGTH-CURPOS ;
        *action = PUSHSNAP ;
        is2++ ;
        LENGTH = remainingLength ;
        CURPOS = 0 ;
        adBinomial_setCut() ;
        (*step)-- ;
      } else {
        *action = ADVANCE ;
        CURPOS++ ;
      }
    }
// printf("NEXT BINOMIAL ACTION:%s\n", actionText[*action]) ;
    return 1 ;
  }
}
      
void adBinomial_resize() {
  int step = 1, i ;
  for (i=BOTTOM+1 ; i<=is2 ; ++i) step += stack2[3*i+1] ;
  printf("Binomial iteration exits on step %i before expected %i\n",
         step-1, stack2[3*BOTTOM]) ;
  stack2[3*BOTTOM] = step-1 ;
  LENGTH = CURPOS ;
  --(CURPOS) ;
}

/****************** INTERFACE CALLED FROM FORTRAN *******************/

void adbinomial_init_(int *length, int *nbSnap, int *firstStep) {
  adBinomial_init(*length, *nbSnap, *firstStep) ;
}

int adbinomial_next_(int *action, int *step) {
  return adBinomial_next(action, step) ;
}

void adbinomial_resize_() {
  adBinomial_resize() ;
}
