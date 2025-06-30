/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2024 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 * 
 * This package gathers functions used to profile a piece of (adjoint differentiated) code in terms of memory and time. 
 */

// Compile with with -D_ADPROFILETRACE=<num> to trace profile computations about the checkpoint
//  whose ckpRank (from 1 up) is <num>, or about all checkpoints if passing -D_ADPROFILETRACE=-1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>

#include <assert.h>

#include "adProfile.h"
#include "adStack.h"

typedef uint64_t stack_size_t;
typedef int64_t  delta_size_t;

typedef clock_t cycle_time_t;
cycle_time_t mytime_() {
        return clock();
}

/** The maximum number of code checkpoint locations that we can handle */
#define MAX_NB_CHECKPOINTS 2000
/** Elementary cost/benefit data of the decision to NOT checkpoint a given static checkpoint location. */
typedef struct _CostBenefit{
  unsigned int ckpRank; // Rank (>=1) of the static checkpoint location this cost/benefit applies to. 
  cycle_time_t DeltaT; // The spared runtime benefit of NOT checkpointing this location.
                       // It is always <=0 (i.e. run-time decreases), and is here stored as unsigned.
  delta_size_t DeltaPk; // The extra peak memory cost of NOT checkpointing this location.
  delta_size_t DeltaTn; // Same as DeltaPk, only about the peak memory at the next turn.
  struct _CostBenefit *next; // Following CostBenefit infos, about other checkpoint locations.
} CostBenefit;

/** One element in a list of runtime sub-checkpoints, immediately below
 * one node/level in the tree of nested runtime checkpoints.
 * The list is ordered from the most recent (downstream) checkpoint till the most upstream one. */
typedef struct _LevelCkpTree {
  unsigned int ckpRank; // The rank (>=1) of the static checkpoint that is now encountered
  cycle_time_t t1t2;     // Time to run primal code+write and read snapshot for this checkpoint.
  cycle_time_t basisCkpT; // Time at some reference point, used to compute time differences.
#ifdef _ADPROFILETRACE
  stack_size_t SnpBefore; // Stack size just before writing the snapshot (only for debug).
#endif
  stack_size_t Snp; // Stack size just after writing the snapshot.
  struct _LevelCkpTree *previous; // Remaining list of runtime checkpoints at this level.
} LevelCkpTree;

/** One node/level in the tree of nested runtime checkpoints */
typedef struct _CkpTree {
  LevelCkpTree *levelCheckpoints; // Ordered list of dynamic sub-checkpoints encountered fwd and not yet bwd, at this level.
  CostBenefit *costBenefits; // Ordered list of cost/benefit information at this level, about each static checkpoint encountered.
  stack_size_t Pk; // Peak stack size reached during reversal of this level.
  stack_size_t Tn; // Peak stack size reached at the next turn.
  LevelCkpTree *repeatLevel; // for management of start/reset/end Repeat
  struct _CkpTree *above; // checkpoints level above (upper un the checkpoints tree)
  struct _CkpTree *below; // checkpoints level below (lower in the checkpoints tree). Kept to avoid malloc/free's.
} CkpTree;

/** The main global tree of active checkpoints at some present time during the profiling run. */
CkpTree initialCkpTree = {.levelCheckpoints=NULL, .costBenefits=NULL,
                            .Pk=0, .Tn=0, .repeatLevel=NULL,
                            .above=NULL, .below=NULL};
CkpTree *curCheckpoints = &initialCkpTree;

/************************* UTILITIES ************************/

/** Describes one (static) checkpoint location in the code source, with a unique ckpRank attached to it.
 * A (yes/no) checkpointing decision applies to one static checkpoint location globally. */
typedef struct {
  char* fname; // Name of the file in which the checkpoint appears
  unsigned int line_nb; // Line number in the original code where the checkpoint appears
  char* callee_name; // Name of the procedure being called (when checkpoint is on a call)
  unsigned int occurrences; // Counter of the number of dynamic occurrences of this static checkpoint.
  unsigned int ckpRank; // Unique integer rank (from 1 up) attached to this checkpoint location.
} CkpLocation;

/** The array of all static checkpoint locations encountered so far */
CkpLocation* allCkpLocations[MAX_NB_CHECKPOINTS] = {};

/** Get or create the unique ckpRank (>=1) associated to the given static checkpoint location,
 * which is the combination of the file name and the line number in this file at which this
 * checkpoint is located. When this checkpoint is at a call site, callname should be passed
 * the procedure name, otherwise NULL. */
unsigned int getCheckpointLocationRank(char* callname, char* filename, unsigned int lineno) {
  unsigned int i = 1;
  int found = 0;
  while (i<MAX_NB_CHECKPOINTS && !found) {
    if (allCkpLocations[i]==NULL) {
      allCkpLocations[i] = (CkpLocation *)malloc(sizeof(CkpLocation));
      allCkpLocations[i]->ckpRank = i;
      allCkpLocations[i]->fname = (char*)malloc((1+strlen(filename))*sizeof(char));
      strcpy(allCkpLocations[i]->fname, filename);
      allCkpLocations[i]->line_nb = lineno;
      if (callname) {
        allCkpLocations[i]->callee_name = (char*)malloc((1+strlen(callname))*sizeof(char));
        strcpy(allCkpLocations[i]->callee_name, callname);
      } else {
        allCkpLocations[i]->callee_name = NULL;
      }
      allCkpLocations[i]->occurrences = 0;
      found = 1;
    } else if (allCkpLocations[i]->line_nb==lineno
               && 0==strcmp(allCkpLocations[i]->fname, filename)
               && 0==strcmp(allCkpLocations[i]->callee_name, callname)) {
      found = 1;
    }
    if (!found) ++i;
  }
  return i;
}

LevelCkpTree* newLevelCkp(unsigned int ckpRank, LevelCkpTree* levelCheckpoints) {
  LevelCkpTree* additionalLevelCkp = (LevelCkpTree*)malloc(sizeof(LevelCkpTree));
  additionalLevelCkp->ckpRank = ckpRank;
  additionalLevelCkp->t1t2 = 0;
  additionalLevelCkp->basisCkpT = 0;
#ifdef _ADPROFILETRACE
  additionalLevelCkp->SnpBefore = 0;
#endif
  additionalLevelCkp->Snp = 0;
  additionalLevelCkp->previous = levelCheckpoints;
  return additionalLevelCkp;
}

void releaseLevelCkp(CkpTree* checkpoints) {
  LevelCkpTree* levelCkps = checkpoints->levelCheckpoints;
  checkpoints->levelCheckpoints = levelCkps->previous;
  if (checkpoints->repeatLevel==NULL) { // Don't free if we are in Repeat mode!
    free(levelCkps);
  }
}

CkpTree* openNewCkp(CkpTree* checkpoints) {
  CkpTree* newCkp = checkpoints->below;
  if (!newCkp) {
    newCkp = (CkpTree*)malloc(sizeof(CkpTree));
    checkpoints->below = newCkp;
    newCkp->repeatLevel = NULL;
    newCkp->above = checkpoints;
    newCkp->below = NULL;
  }
  newCkp->levelCheckpoints = NULL;
  newCkp->costBenefits = NULL;
  newCkp->Pk = 0;
  newCkp->Tn = 0;
  return newCkp;
}

CkpTree* closeCkp(CkpTree* checkpoints) {
  return checkpoints->above;
}

/** Advances in the given (address of a) chain of CostBenefit's, till finding a
 * CostBenefit about the checkpoint indexed as "ckpRank". This chain is ordered by growing ckpRank's.
 * If no such CostBenefit is found, creates and inserts a new, empty one, at the correct location.
 * This insertion modifies the given chain of CostBenefit's by side-effect.
 * Returns the address of the sub-chain of CostBenefit that starts with the one found or created. */
CostBenefit** getSetCostBenefit(CostBenefit **toCostBenefits, unsigned int ckpRank) {
  while (*toCostBenefits!=NULL && (*toCostBenefits)->ckpRank < ckpRank) {
    toCostBenefits = &((*toCostBenefits)->next);
  }
  if (*toCostBenefits==NULL || (*toCostBenefits)->ckpRank > ckpRank) {
    CostBenefit* additionalCostBenefits = (CostBenefit*)malloc(sizeof(CostBenefit));
    additionalCostBenefits->ckpRank = ckpRank;
    additionalCostBenefits->DeltaT = 0;
    additionalCostBenefits->DeltaPk = 0;
    additionalCostBenefits->DeltaTn = 0;
    additionalCostBenefits->next = *toCostBenefits;
    *toCostBenefits = additionalCostBenefits;
  }
  return toCostBenefits;
}

#ifdef _ADPROFILETRACE

int depthCkp(CkpTree* checkpoints) {
  int i = 0;
  while (checkpoints) {
    ++i;
    checkpoints = checkpoints->above;
  }
  return i;
}

int costBenefitsLength(CostBenefit* costBenefits) {
  int i = 0;
  while (costBenefits) {
    ++i;
    costBenefits = costBenefits->next;
  }
  return i;
}

void dumpCostBenefits(CostBenefit *costBenefits) {
  while (costBenefits) {
    printf(" %u:(D_t:%"PRIu64"; D_Tn:%"PRId64"; D_Pk:%"PRId64")",
           costBenefits->ckpRank, costBenefits->DeltaT, costBenefits->DeltaTn, costBenefits->DeltaPk);
    costBenefits = costBenefits->next;
  }
}

#endif

/****************** MAIN PUBLISHED PRIMITIVES ****************/

void adProfileAdj_SNPWrite(char* callname, char* filename, unsigned int lineno) {
  unsigned int ckpRank = getCheckpointLocationRank(callname, filename, lineno);
  curCheckpoints->levelCheckpoints = newLevelCkp(ckpRank, curCheckpoints->levelCheckpoints);
#ifdef _ADPROFILETRACE
  curCheckpoints->levelCheckpoints->SnpBefore = adStack_getCurrentStackSize();
#endif
  curCheckpoints->levelCheckpoints->basisCkpT = mytime_();
}

void adProfileAdj_beginAdvance(char* callname, char* filename, unsigned int lineno) {
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  curCheckpoints->levelCheckpoints->Snp = adStack_getCurrentStackSize();
}

void adProfileAdj_endAdvance(char* callname, char* filename, unsigned int lineno) {
  curCheckpoints->levelCheckpoints->t1t2 = mytime_() - curCheckpoints->levelCheckpoints->basisCkpT;
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
}

void adProfileAdj_SNPRead(char* callname, char* filename, unsigned int lineno) {
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  curCheckpoints->levelCheckpoints->basisCkpT = mytime_();
}

void adProfileAdj_beginReverse(char* callname, char* filename, unsigned int lineno) {
  curCheckpoints->levelCheckpoints->t1t2 += mytime_() - curCheckpoints->levelCheckpoints->basisCkpT;
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  ++(allCkpLocations[curCheckpoints->levelCheckpoints->ckpRank]->occurrences);
  curCheckpoints = openNewCkp(curCheckpoints);
}

void adProfileAdj_endReverse(char* callname, char* filename, unsigned int lineno) {
  assert(curCheckpoints->levelCheckpoints==NULL);
  CostBenefit *costBenefits_C = curCheckpoints->costBenefits;
  stack_size_t Pk_C = curCheckpoints->Pk;
  stack_size_t Tn_C = curCheckpoints->Tn;
  // Special case if "callname" is an external, "callname_B" did not call adProfileAdj_turn,
  //  and curCheckpoints->Pk is still 0 ! Set it to the current stack size:
  if (Pk_C==0) Pk_C = adStack_getCurrentStackSize();
  curCheckpoints = closeCkp(curCheckpoints);
  // Now merge the additional CostBenefit's from the checkpointed fragment C into the current CostBenefit's:
  unsigned int rankOfC = curCheckpoints->levelCheckpoints->ckpRank;
  stack_size_t Pk_D = curCheckpoints->Pk;
  stack_size_t Tn_D = curCheckpoints->Tn;
  int mergedC = 0;
  unsigned int rankOfX;
  cycle_time_t DeltaT_CX;
  delta_size_t DeltaPk_CX;
  delta_size_t DeltaTn_CX;
#ifdef _ADPROFILETRACE
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==rankOfC) {
    printf("ENDREVERSE OF %s::%i CALL %s, DEPTH:%i, SIZEBEFORESNP:%"PRId64", SIZEAFTERSNP:%"PRId64"\n",
           allCkpLocations[rankOfC]->fname,
           allCkpLocations[rankOfC]->line_nb,
           (allCkpLocations[rankOfC]->callee_name ? allCkpLocations[rankOfC]->callee_name : "MANUAL"),
           depthCkp(curCheckpoints),
           curCheckpoints->levelCheckpoints->SnpBefore,
           curCheckpoints->levelCheckpoints->Snp);
    printf("  BEFORE (Tn_D:%"PRId64" bytes, Pk_D:%"PRId64" bytes) :", Tn_D, Pk_D);
    dumpCostBenefits(curCheckpoints->costBenefits);
    printf("\n  ADDING (Tn_C:%"PRId64" bytes, Pk_C:%"PRId64" bytes) :", Tn_C, Pk_C);
    dumpCostBenefits(costBenefits_C);
    printf("\n");
  }
  int nbProfilesBefore = costBenefitsLength(curCheckpoints->costBenefits);
  int nbProfilesAdded = costBenefitsLength(costBenefits_C);
#endif
  CostBenefit **toCostBenefits_D = &(curCheckpoints->costBenefits);
  // compute merged peak size (Tn_D does not change):
  if (Pk_C>Pk_D) curCheckpoints->Pk = Pk_C;
  // merge C's list of costs/benefits into D's:
  while (costBenefits_C || !mergedC) {
    rankOfX = (costBenefits_C ? costBenefits_C->ckpRank : 0);
#ifdef _ADPROFILETRACE
    if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==rankOfC) {
      printf("    mergedC:%i next additional:%u,t:%"PRIu64",Tn%"PRId64",Pk:%"PRId64"\n",
             mergedC, rankOfX,
             (costBenefits_C?costBenefits_C->DeltaT:0),
             (costBenefits_C?costBenefits_C->DeltaTn:0),
             (costBenefits_C?costBenefits_C->DeltaPk:0));
  }
#endif
    if (!mergedC
        && (costBenefits_C==NULL || rankOfC < rankOfX)) {
      // it is time to incorporate costBenefit of C alone, even if it not the candidate for merge X, to preserve ordering:
      rankOfX = rankOfC;
      DeltaT_CX = curCheckpoints->levelCheckpoints->t1t2;
      DeltaPk_CX = 0;
      DeltaTn_CX = 0;
      mergedC = 1;
    } else if (rankOfX == rankOfC) {
      // it is time to incorporate constBenefit of C, together with candidate for merge X (which is the same):
      DeltaT_CX = costBenefits_C->DeltaT + curCheckpoints->levelCheckpoints->t1t2;
      DeltaPk_CX = costBenefits_C->DeltaPk ;
      DeltaTn_CX = costBenefits_C->DeltaTn;
      CostBenefit* tmp = costBenefits_C->next;
      free(costBenefits_C);
      costBenefits_C = tmp;
      mergedC = 1;
    } else {
      // else incorporate the candidate for merge X, alone:
      DeltaT_CX = costBenefits_C->DeltaT;
      DeltaPk_CX = costBenefits_C->DeltaPk;
      DeltaTn_CX = costBenefits_C->DeltaTn;
      CostBenefit* tmp = costBenefits_C->next;
      free(costBenefits_C);
      costBenefits_C = tmp;
    }
    toCostBenefits_D = getSetCostBenefit(toCostBenefits_D, rankOfX);
#ifdef _ADPROFILETRACE
    if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==rankOfC) {
      printf("    X:%u DeltaPk:%"PRId64" DeltaPk_C:%"PRId64"\n", rankOfX, (*toCostBenefits_D)->DeltaPk, DeltaPk_CX);
    }
#endif
    delta_size_t offsetOfNoCkp_C =
        (rankOfX==rankOfC ? Tn_C + DeltaTn_CX - curCheckpoints->levelCheckpoints->Snp : 0);
#ifdef _ADPROFILETRACE
    if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==rankOfC) {
      printf("      %"PRIu64"  %"PRId64" %"PRIu64"\n", Tn_C, DeltaTn_CX, curCheckpoints->levelCheckpoints->Snp) ;
      printf("      %"PRId64" += %"PRId64"\n", (*toCostBenefits_D)->DeltaTn, offsetOfNoCkp_C);
    }
#endif
    (*toCostBenefits_D)->DeltaTn += offsetOfNoCkp_C;
    stack_size_t newPk_D = Pk_D + (*toCostBenefits_D)->DeltaPk + offsetOfNoCkp_C;
    stack_size_t newPk_C = Pk_C + DeltaPk_CX;
    (*toCostBenefits_D)->DeltaPk = (newPk_D>newPk_C ? newPk_D : newPk_C) - curCheckpoints->Pk;
    (*toCostBenefits_D)->DeltaT += DeltaT_CX;
    toCostBenefits_D = &((*toCostBenefits_D)->next);
  }
#ifdef _ADPROFILETRACE
  int nbProfilesAfter = costBenefitsLength(curCheckpoints->costBenefits);
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==rankOfC) {
    printf("  ==> MERGE %i COSTBENEFITS INTO %i --> %i\n",
           nbProfilesAdded, nbProfilesBefore, nbProfilesAfter);
    printf("  >AFTER (Tn_CD:%"PRId64" bytes, Pk_CD:%"PRId64" bytes) :",
           curCheckpoints->Tn, curCheckpoints->Pk);
    dumpCostBenefits(curCheckpoints->costBenefits);
    printf("\n");
  }
#endif
  // merge finished.
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  releaseLevelCkp(curCheckpoints);
}

void adProfileAdj_turn(char* callname, char* filename) {
  stack_size_t curPeak = adStack_getCurrentStackSize();
  if (curPeak > curCheckpoints->Pk) {
    curCheckpoints->Pk = curPeak;
    curCheckpoints->Tn = curPeak;
  }
}

void adProfileAdj_startRepeat() {
  curCheckpoints->repeatLevel = curCheckpoints->levelCheckpoints;
}

void adProfileAdj_resetRepeat() {
  curCheckpoints->levelCheckpoints = curCheckpoints->repeatLevel;
}

void adProfileAdj_endRepeat() {
  while (curCheckpoints->repeatLevel != curCheckpoints->levelCheckpoints) {
    // Now that we are no longer in Repeat mode, free preserved LevelCkpTree's:
    LevelCkpTree *previousElem = curCheckpoints->repeatLevel->previous;
    free(curCheckpoints->repeatLevel);
    curCheckpoints->repeatLevel = previousElem;
  }
  curCheckpoints->repeatLevel = NULL;
}

typedef struct _SortedCostBenefit {
  unsigned int ckpRank; 
  cycle_time_t DeltaT;
  delta_size_t DeltaPk;
  float ratio;
  struct _SortedCostBenefit* next;
} SortedCostBenefit;

int callNameComesAfter(unsigned int ckpRank1, unsigned int ckpRank2) {
  int comparison = strcmp(allCkpLocations[ckpRank1]->callee_name, allCkpLocations[ckpRank2]->callee_name);
  return (comparison>0 || (comparison==0 && allCkpLocations[ckpRank1]->line_nb > allCkpLocations[ckpRank2]->line_nb));
}

void showOneCostBenefit(SortedCostBenefit* sortedCostBenefits, FILE* fp) {
  CkpLocation* ckpLocation = allCkpLocations[sortedCostBenefits->ckpRank];
  cycle_time_t DeltaT = (sortedCostBenefits->DeltaT / CLOCKS_PER_SEC)*1000; // DeltaT in milliSeconds.
  int DeltaTsec = DeltaT/1000;
  int DeltaTmillisec = DeltaT%1000;
  fprintf(fp, "  - Time gain -%2d.%03d s.", DeltaTsec, DeltaTmillisec);
  delta_size_t DeltaPk = sortedCostBenefits->DeltaPk;
  if (DeltaPk<0) {
    fprintf(fp, " and peak memory gain %"PRId64"b", DeltaPk);
  } else if (DeltaPk==0) {
    fprintf(fp, " at peak memory cost zero");
  } else {
    int DeltaPkMb = DeltaPk/1000000;
    int DeltaPkb = DeltaPk%1000000;
    fprintf(fp, " at peak memory cost %4d.%06d Mb", DeltaPkMb, DeltaPkb);
  }
  if (ckpLocation->callee_name) {
    fprintf(fp, " for call %s (%u times), at", ckpLocation->callee_name, ckpLocation->occurrences);
  } else {
    fprintf(fp, " for checkpoint (%u times) starting at", ckpLocation->occurrences);
  }
  fprintf(fp, " location#%u: line %u of file %s\n", sortedCostBenefits->ckpRank, ckpLocation->line_nb, ckpLocation->fname);
  sortedCostBenefits = sortedCostBenefits->next;
}

void adProfileAdj_showProfiles() {
    adProfileAdj_showProfilesFile(NULL);
}

void adProfileAdj_showProfilesFile(char* filename) {
  SortedCostBenefit* sortedCostBenefitsN = NULL;
  SortedCostBenefit* sortedCostBenefitsZ = NULL;
  SortedCostBenefit* sortedCostBenefitsP = NULL;
  CostBenefit *costBenefits = initialCkpTree.costBenefits;
  SortedCostBenefit** toSortedCostBenefits;
  float ratio;
  FILE* fpout = stdout;
  if( filename!=NULL && strlen(filename)>0 ) {
      fpout = fopen(filename, "w");
      if( fpout==NULL ) {
          fprintf(stderr, "ERROR::***%s*** was not opened properly", filename);
          exit(1);
      }
  }
  while (costBenefits) {
    unsigned int newCkpRank = costBenefits->ckpRank;
    if (costBenefits->DeltaPk <= 0) {
      // Sort negative memory costs and zero memory costs in two lists ordered by routine name then location rank:
      toSortedCostBenefits = (costBenefits->DeltaPk==0 ? &sortedCostBenefitsZ : &sortedCostBenefitsN);
      ratio = 0.0;
      while (*toSortedCostBenefits && callNameComesAfter(newCkpRank, (*toSortedCostBenefits)->ckpRank)) {
        toSortedCostBenefits = &((*toSortedCostBenefits)->next);
      }
    } else {
      // Sort positive memory costs by decreasing time/memory benefit:
      toSortedCostBenefits = &sortedCostBenefitsP;
      ratio = ((float)costBenefits->DeltaT) / ((float)costBenefits->DeltaPk);
      while (*toSortedCostBenefits && ratio<(*toSortedCostBenefits)->ratio) {
        toSortedCostBenefits = &((*toSortedCostBenefits)->next);
      }
    }
    SortedCostBenefit *newSorted = (SortedCostBenefit*)malloc(sizeof(SortedCostBenefit));
    newSorted->ckpRank = newCkpRank;
    newSorted->DeltaT = costBenefits->DeltaT;
    newSorted->DeltaPk = costBenefits->DeltaPk;
    newSorted->ratio = ratio;
    newSorted->next = *toSortedCostBenefits;
    *toSortedCostBenefits = newSorted;
    costBenefits = costBenefits->next;
  }
  fprintf(fpout, "PEAK STACK:%"PRId64" bytes\n", initialCkpTree.Pk);
  fprintf(fpout, "SUGGESTED NOCHECKPOINTs:\n");
  fprintf(fpout, " * Peak memory gain:\n");
  while (sortedCostBenefitsN) {
      showOneCostBenefit(sortedCostBenefitsN, fpout);
    sortedCostBenefitsN = sortedCostBenefitsN->next;
  }
  fprintf(fpout, " * Peak memory neutral:\n");
  while (sortedCostBenefitsZ) {
      showOneCostBenefit(sortedCostBenefitsZ, fpout);
    sortedCostBenefitsZ = sortedCostBenefitsZ->next;
  }
  fprintf(fpout, " * Peak memory cost:\n");
  while (sortedCostBenefitsP) {
      showOneCostBenefit(sortedCostBenefitsP, fpout);
    sortedCostBenefitsP = sortedCostBenefitsP->next;
  }
  if( fpout!=stdout ) {
      fclose(fpout);
      printf("profiling information written to ***%s***\n\n", filename);
  }
}

/****************** INTERFACE CALLED FROM FORTRAN *******************/

void adprofileadj_snpwrite_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_SNPWrite(callname, filename, *lineno);
}

void adprofileadj_beginadvance_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_beginAdvance(callname, filename, *lineno);
}

void adprofileadj_endadvance_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_endAdvance(callname, filename, *lineno);
}

void adprofileadj_snpread_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_SNPRead(callname, filename, *lineno);
}

void adprofileadj_beginreverse_(char* callname, char* filename, unsigned int *lineno){
  adProfileAdj_beginReverse(callname, filename, *lineno);
}

void adprofileadj_endreverse_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_endReverse(callname, filename, *lineno);
}

void adprofileadj_turn_(char* callname, char* filename) {
  adProfileAdj_turn(callname, filename);
}

void adprofileadj_startrepeat_() {
  adProfileAdj_startRepeat();
}

void adprofileadj_resetrepeat_() {
  adProfileAdj_resetRepeat();
}

void adprofileadj_endrepeat_() {
  adProfileAdj_endRepeat();
}

void adprofileadj_showprofilesfile_(char* filename) {
    char* ftnfile;
    int len;
    len = strlen(filename) + 1;
    ftnfile = (char*) malloc( len*sizeof(char) );
    strcpy(ftnfile, filename);
    ftnfile[len-1] = '\0';
    adProfileAdj_showProfilesFile(ftnfile);
    free( ftnfile );
}

void adprofileadj_showprofiles_() {
  adProfileAdj_showProfiles();
}
