#ifndef ADPROFILE_LOADED
#define ADPROFILE_LOADED

/** To be called just before writing a snapshot */
void adProfileAdj_SNPWrite(char* callname, char* filename, unsigned int lineno);
/** To be called just after writing a snapshot and before the duplicated primal call */
void adProfileAdj_beginAdvance(char* callname, char* filename, unsigned int lineno);
/** To be called just after the duplicated primal call */
void adProfileAdj_endAdvance(char* callname, char* filename, unsigned int lineno);

/** To be called just before reading a snapshot */
void adProfileAdj_SNPRead(char* callname, char* filename, unsigned int lineno);
/** To be called just after reading a snapshot and before the checkpointed adjoint call */
void adProfileAdj_beginReverse(char* callname, char* filename, unsigned int lineno);
/** To be called just after the checkpointed adjoint call */
void adProfileAdj_endReverse(char* callname, char* filename, unsigned int lineno);

/** To be called at each turn point, when stack size reaches a peak */
void adProfileAdj_turn(char* callname, char* filename);

/** To be called along with adStack_startRepeat() in FixedPoint adjoint code */
void adProfileAdj_startRepeat() ;
/** To be called along with adStack_resetRepeat() in FixedPoint adjoint code */
void adProfileAdj_resetRepeat() ;
/** To be called along with adStack_endRepeat() in FixedPoint adjoint code */
void adProfileAdj_endRepeat() ;

/** Displays the estimated costs/benefits. To be called at the end of the profiled execution. */
void adProfileAdj_showProfiles() ;
void adProfileAdj_showProfilesFile(char* filename) ;


#endif // ADPROFILE_LOADED
