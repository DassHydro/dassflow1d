/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2021 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "adContext.h"
#include "adComplex.h"
#include "json.h"


typedef enum dbad_mode  { TANGENT=1, ADJOINT=2 } dbad_mode_t;
typedef enum seed_value { SEED_DEFAULT=0, SEED_OFF=-1 } seed_value_t;
typedef enum var_type   { VAR_SCALAR=1, VAR_ARRAY=2 } var_type_t;
typedef struct seed_spec_t { int ndir; int* dir_lst; } seed_spec_t;
static int dbad_mode, dbad_phase ;
static double dbad_ddeps = 1.e-6 ;
static double dbad_condensed_val, dbad_condensed_tgt, dbad_condensed_adj ;

// 
static double dbad_seed = 0.137 ;
static double dbad_currentSeed;
// settings to be loaded from JSON file (if present in workdir)
// { 
//     "eps":1.e-8, 
//     "seed":0.137, 
//     "only_indep":[
//                   {"x":[1,8]}, 
// 	             "y"
//                  ], 
//     "only_dep": [
//                  "z",
//                  {"z1":[3]}
//                 ] 
// } 

// default values for
// - perturbation (eps) for divided differences
// - seeds for independents/dependents
// may be changed by user specified json file.
const static char* jsonfile = "dbad.json";

typed(json_object) * glob_json_obj = NULL; 
typed(json_element)  glob_json_element;


const int extdebug = 0; //extended debug output (potentially helpful in case of pitfalls)


/***************************************************/
/*        SEED handling
/**************************************************/
static seed_spec_t* seed_spec_alloc(int ndir) {
    seed_spec_t* pseed = malloc( sizeof(seed_spec_t) );
    pseed->ndir = ndir;
    if(ndir>0) {
	pseed->dir_lst = malloc(ndir*sizeof(int));
    }
    else {
	pseed->dir_lst = NULL;
    }
}

static void seed_spec_dispose(seed_spec_t* pseed) {
    if( pseed!=NULL ) {
	if(pseed->dir_lst!=NULL) free(pseed->dir_lst);
	free(pseed);
    }
}

static void seed_spec_dump(seed_spec_t* pseed) {
    if( pseed==NULL ) {
        printf("(nil)");
    }
    else if( pseed->ndir==SEED_DEFAULT ) {
        printf("(default)");
    }
    else if( pseed->ndir==SEED_OFF) {
        printf("(off)");
    }
    else {
        printf("(indices: ");
        for(int i=0; i<pseed->ndir; ++i) {
            printf("%d, ", pseed->dir_lst[i]);
        }
        printf(")");
    }
}


/***************************************************/
/*        reading file into char buffer
/* Note:
/* caller is reponsible to release allocated buffer
/*
/**************************************************/
static const char *read_file(const char *path) {

  FILE *file = fopen(path, "r");

  if (file == NULL) {
      //fprintf(stderr, "Expected file \"%s\" not found", path);
      return NULL;
  }

  fseek(file, 0, SEEK_END);
  long len = ftell(file);
  fseek(file, 0, SEEK_SET);
  char *buffer = malloc(len + 1);

  if (buffer == NULL) {
    fprintf(stderr, "Unable to allocate memory for file");
    fclose(file);
    return NULL;
  }

  fread(buffer, 1, len, file);
  buffer[len] = '\0';

  fclose(file);

  return (const char *)buffer;
}

/***************************************************/
/*        JSON file reading
/**************************************************/
static const char* json_type_str(typed(json_element_type) type) {
    switch(type) {
        case JSON_ELEMENT_TYPE_NULL:
            return "null";
            break;
        case JSON_ELEMENT_TYPE_BOOLEAN:
            return "bool";
            break;
        case JSON_ELEMENT_TYPE_STRING:
            return "string";
            break;
        case JSON_ELEMENT_TYPE_NUMBER:
            return "number";
            break;
        case JSON_ELEMENT_TYPE_OBJECT:
            return "object";
            break;
        case JSON_ELEMENT_TYPE_ARRAY:
            return "array";
            break;
        default:
            fprintf(stderr, "unexpected json_element_type");
            return NULL;
            break;
    }
}


static void jsonelem_get_double(typed(json_element) elem, double* value, int* ok) {
    *ok = 0;
    if( elem.type!=JSON_ELEMENT_TYPE_NUMBER ) {
        fprintf(stderr, "unexpected type -->%s<-- when reading double\n", json_type_str(elem.type));
        return;
    }
    else {
        if( elem.value.as_number.type==JSON_NUMBER_TYPE_DOUBLE ) {
            *value = elem.value.as_number.value.as_double;
            *ok = 1;
        }
        else {//JSON_NUMBER_TYPE_LONG
            *value = (double) elem.value.as_number.value.as_long;
            *ok = 1;
        }
    }
}


static void jsonelem_get_i4(typed(json_element) elem, int* value, int* ok) {
    *ok = 0;
    if( elem.type!=JSON_ELEMENT_TYPE_NUMBER ) {
        fprintf(stderr, "unexpected type -->%s<-- when reading double\n", json_type_str(elem.type));
        return;
    }
    else {
        if( elem.value.as_number.type==JSON_NUMBER_TYPE_DOUBLE ) {
	    fprintf(stderr, "json double not expected here\n");
	    exit(1);
        }
        else {//JSON_NUMBER_TYPE_LONG
	    typed(json_number_long) lvalue = elem.value.as_number.value.as_long;
            *value = (int) lvalue; //TODO::check that lvalue does not exceed range of integers
            *ok = 1;
        }
    }
}

static typed(json_object) * jsonobj_find_obj(typed(json_object) *obj, const char* key) {
    typed(json_object) * retobj = NULL;
    if (dbad_phase>=999) {
        printf("%s:: key -->%s<--\n", __func__, key);
    }
    result(json_element) key_result = json_object_find(obj, key);
    if( result_is_err(json_element)(&key_result)) {
	// key not found, nothing to change
	return  retobj;
    }
    typed(json_element) elem = result_unwrap(json_element)(&key_result);
    if( elem.type!=JSON_ELEMENT_TYPE_OBJECT ) {
	fprintf(stderr, "unexpected type -->%s<-- for KEY -->%s<--\n",
                json_type_str(elem.type), key);
    }
    else {
	retobj = elem.value.as_object;
    }

    return retobj;
}


static void jsonobj_read_double(typed(json_object) *obj, const char* key, double* value, int* ok) {
    *ok = 0;
    if(dbad_phase>=999) {
        printf("%s:: json_obj=%p\n",__func__, obj);
    }
    result(json_element) key_result = json_object_find(obj, key);
    if( ! result_is_err(json_element)(&key_result)) {
        typed(json_element) elem = result_unwrap(json_element)(&key_result);
        jsonelem_get_double(elem, value, ok);
    }
    else {
        // key not found or other error
    }
}


static double* jsonobj_read_double_lst(typed(json_object) *obj, const char* key, int* cnt) {
    *cnt = 0;
    double* value_lst = NULL;
    double cur_value;
    int ok;

    result(json_element) key_result = json_object_find(obj, key);

    if( result_is_err(json_element)(&key_result)) {
        return value_lst;
    }
    //
    typed(json_element) elem = result_unwrap(json_element)(&key_result);
    if( elem.type!=JSON_ELEMENT_TYPE_ARRAY ) {
        fprintf(stderr, "unexpected type -->%s<-- when reading value list\n",
                json_type_str(elem.type));
    }
    else {
        *cnt = elem.value.as_array->count;
        value_lst = malloc(*cnt*sizeof(double));
        for(int j = 0; j <*cnt; ++j) {
            typed(json_element) cur_elem = elem.value.as_array->elements[j];
            jsonelem_get_double(cur_elem, &cur_value, &ok);
            if( ok ) {
                value_lst[j] = cur_value;
            }
            else {
                //todo message
                free( value_lst );
                value_lst = NULL;
                break;
            }
        }
    }

    return value_lst;
}

static seed_spec_t *
jsonobj_read_seed(typed(json_object) *obj, const char* key, const var_type_t vtyp) {
    seed_spec_t * seed_spec = NULL;
    int cnt = 0;
    int cur_value;
    int ok;

    // lookup the 'key'
    result(json_element) key_result = json_object_find(obj, key);

    typed(json_element) elem;
    if( result_is_err(json_element)(&key_result)) { // key not found
        seed_spec = seed_spec_alloc(SEED_OFF);
        return seed_spec;// key not found
    }
    else {
        typed(json_element) elem = result_unwrap(json_element)(&key_result);
    }
    if(extdebug) printf("EXTDEBUG: KEY=%s  type -->%s<--\n", key, json_type_str(elem.type));
    if( elem.type==JSON_ELEMENT_TYPE_BOOLEAN ) {
        typed(json_boolean) use_flag = elem.value.as_boolean;
        if( use_flag ) {
            seed_spec = seed_spec_alloc(SEED_DEFAULT);
        }
        else {//should not be used
            if( dbad_phase>=99) {
                printf("disabled directional derivative for variable -->%s<--\n", key);
            }
            seed_spec = seed_spec_alloc(SEED_OFF);
            return seed_spec;
        }
    }
    else if( elem.type==JSON_ELEMENT_TYPE_ARRAY ) {
        if( vtyp!=VAR_ARRAY ) {
            fprintf(stderr, "ERROR::array like specification for non array variable -->%s<--\n",
                    key);
            seed_spec = seed_spec_alloc(SEED_DEFAULT);
        }
        else {
            cnt = elem.value.as_array->count;
            seed_spec = seed_spec_alloc(cnt);
            if(extdebug) printf("EXTDEBUG::detected ndir=%d for -->%s<--\n", cnt, key);
            for(int j = 0; j <cnt; ++j) {
                typed(json_element) cur_elem = elem.value.as_array->elements[j];
                jsonelem_get_i4(cur_elem, &cur_value, &ok);
                if( ok ) {
                    seed_spec->dir_lst[j] = cur_value;
                }
                else {
                    fprintf(stderr, "ERROR when reading seed indices for KEY -->%s<--\n",
                            key);
                    seed_spec_dispose(seed_spec);
                    seed_spec = seed_spec_alloc(SEED_DEFAULT);
                    break;
                }
            }
        }
    }
    else {
        fprintf(stderr, "unexpected type -->%s<-- when reading KEY -->%s<--\n",
                json_type_str(elem.type), key);
        exit(EXIT_FAILURE);
    }

    return seed_spec;
}

static void dbad_json_load() {
    // epsilon/seed potentially may be provided via json file,
    // which takes precedence if present
    const char* json_str = read_file(jsonfile);
    double value;
    double* value_lst;
    int rc;
    if( json_str!=NULL ) {
	// still kept
	result(json_element) element_result = json_parse(json_str);
        free((void *)json_str); // dispose character buffer

        if (result_is_err(json_element)(&element_result)) {
            typed(json_error) error = result_unwrap_err(json_element)(&element_result);
            fprintf(stderr, "Error parsing JSON: %s\n", json_error_to_string(error));
            return;
        }
        else {
            if( dbad_phase>=99 ) {
                printf("...start reading settings from file ***%s***\n", jsonfile);
            }

            // Extract the data (go to global variables)
            glob_json_element = result_unwrap(json_element)(&element_result);
            glob_json_obj     = glob_json_element.value.as_object;
        }
    }
}//dbad_json_load


static void dbad_json_set_ddeps() {
    int ok;
    double json_eps;
    if( glob_json_obj!=NULL ) {
        if( dbad_phase>=9999 ) {
            printf("***%s*** --> calling jsonobj_read_double...\n", __func__);
        }
	jsonobj_read_double(glob_json_obj, "eps", &json_eps, &ok);
        if( dbad_phase>=999 ) {
            if(ok) {
                printf("%s:: returns ok=%d eps=%10.4e\n", __func__, ok, json_eps);
            }
            else {
                printf("%s:: reading eps failed ok=%d\n", __func__, ok);
            }
        }
	if( ok ) {
	    dbad_ddeps = json_eps;
            fprintf(stdout,
                    "-->eps<-- provided by json file.\n");
	}
	else {
            fprintf(stdout,
                    "-->eps<-- not provided by json file, leaving default unchanged.\n");
	}
    }
}

static void dbad_json_set_seed() {
    int ok;
    double json_seed;
    if( glob_json_obj!=NULL ) {
        if( dbad_phase>=9999 ) {
            printf("***%s*** --> calling jsonobj_read_double...\n", __func__);
        }
	jsonobj_read_double(glob_json_obj, "seed", &json_seed, &ok);
        if( dbad_phase>=999 ) {
            if(ok) {
                printf("%s:: returns ok=%d seed=%24.16e\n", __func__, ok, json_seed);
            }
            else {
                printf("%s:: reading seed failed ok=%d\n", __func__, ok);
            }
        }
	if( ok ) {
	    dbad_seed = json_seed;
            fprintf(stdout,
                    "-->seed<-- provided by json file.\n");
	}
	else {
            fprintf(stdout,
                    "-->seed<-- not provided by json file, leaving default unchanged.\n");
	}
    }
}

static seed_spec_t * dbad_list_find_seedspec(typed(json_array) * list,
                                            const char* vname, const var_type_t vtyp) {

    seed_spec_t * seedspec = NULL;
    if( dbad_phase>=999 ) {
        printf("***%s*** for vname -->%s<-- len(list)=%d\n", __func__, vname, list->count);
    }
    for(int i=0; i<list->count; ++i) {
        typed(json_element) elem = list->elements[i]; //current element
        if( elem.type==JSON_ELEMENT_TYPE_STRING ) {
            if( dbad_phase>=999 ) {
                printf("***%s*** i=%d detected var -->%s<-- cmp=%d\n",
                       __func__, i, elem.value.as_string, strcmp(elem.value.as_string,vname));
            }
            if( strcmp(elem.value.as_string,vname)==0 ) {
                seedspec = seed_spec_alloc(SEED_DEFAULT);
                break; //vname found
            }
            else {
                continue; // not found, continue searc
            }
        }
        else if( elem.type==JSON_ELEMENT_TYPE_OBJECT ) { //must provide *single* array like spec.
            typed(json_object) * obj = elem.value.as_object;
            // looking up vname
            result(json_element) vname_result = json_object_find(obj, vname);
            if( result_is_err(json_element)(&vname_result)) {
                continue; // vname not found
            }
            else {
                typed(json_element) elem = result_unwrap(json_element)(&vname_result);
                if( elem.type!=JSON_ELEMENT_TYPE_ARRAY ) {
                    fprintf(stderr, "***%s*** ERROR seed specification for variable -->%s<--"
                            " array type expected but found type ==%s==\n",
                            __func__, vname, json_type_str(elem.type));
                    exit(EXIT_FAILURE);
                }
                else {
                    if( vtyp!=VAR_ARRAY ) {
                        fprintf(stderr, "***%s*** ERRROR array like seed specification for non-array variable -->%s<--, falling back to default!\n", __func__, vname);
                        seedspec = seed_spec_alloc(SEED_DEFAULT);
                        break;
                    }
                    else {
                        typed(json_array) * arr = elem.value.as_array;
                        seedspec = seed_spec_alloc(arr->count);
                        for(int j=0; j<arr->count; ++j) {
                            typed(json_element) cur_elem = elem.value.as_array->elements[j];
                            int cur_value, ok;
                            jsonelem_get_i4(cur_elem, &cur_value, &ok);
                            if(ok) {
                                seedspec->dir_lst[j] = cur_value;
                            }
                            else {
                                fprintf(stderr, "***%s*** ERROR inconsistent seed index specification for variable -->%s<--\n",
                                        __func__, vname);
                            }
                        }
                        break;
                    }
                }
            }
        }
        else {
            fprintf(stderr, "***%s*** ERROR seed specification for variable -->%s<--"
                    " with inconsistent type ==%s==\n", __func__, vname, json_type_str(elem.type));
            exit(EXIT_FAILURE);
        }
    }

    if( seedspec==NULL ) { //'varname' was not found
        seedspec = seed_spec_alloc(SEED_OFF);
    }
    return seedspec;
}


static seed_spec_t * dbad_json_get_seedspec(const char* key, const char* varname, const var_type_t vtyp) {
    // (1) dbad.json missing:                 seed default
    // (2) key not found in dbad.json:        seed default
    // (3) varname not found in key list:     seed off
    // (4) varname found in key list:         seed default (potentially only for selected array components)
    seed_spec_t * seedspec = NULL;

    if(dbad_phase>=999) {
        printf("***%s*** key -->%s<-- varname -->%s<--\n",__func__, key, varname);
    }

    if( glob_json_obj!=NULL ) {
        if(dbad_phase>=999) {
            printf("%s:: glob_json_obj=%p\n",__func__, glob_json_obj);
        }
        // looking up only_indep/only_dep
        result(json_element) key_result = json_object_find(glob_json_obj, key);
        typed(json_element) elem;
        if( result_is_err(json_element)(&key_result)) { // key  not found
            if(dbad_phase>=999) {
                printf("***%s*** key -->%s<-- NOT found\n",__func__, key);
            }
            seedspec = seed_spec_alloc(SEED_DEFAULT);
        }
        else {
            elem = result_unwrap(json_element)(&key_result);
            if( elem.type!=JSON_ELEMENT_TYPE_ARRAY ) {
                fprintf(stderr, "***%s*** ERROR:: -->%s<-- is NOT list type!\n", __func__, key);
                exit(EXIT_FAILURE);
            }
            else {
                seedspec = dbad_list_find_seedspec(elem.value.as_array, varname, vtyp);
            }
        }
    }

    if( dbad_phase>=999 ) {
        printf("%s seed for -->%s<-- ", __func__, varname);
        seed_spec_dump(seedspec);
        printf("\n");
    }

    return seedspec;
}

static seed_spec_t * dbad_json_get_indep(const char* varname, const var_type_t vtyp) {
    // (1) dbad.json missing:                 seed default
    // (2) no 'only_indep' in dbad.json:      seed default
    // (3) varname not found in 'only_indep': seed off
    // (4) varname found in 'only_indep':     seed default (potentially only for selected array components)
    return dbad_json_get_seedspec("only_indep", varname, vtyp);
}

static seed_spec_t * dbad_json_get_dep(const char* varname, const var_type_t vtyp) {
    // (1) dbad.json missing:               seed default
    // (2) no 'only_dep' in dbad.json:      seed default
    // (3) varname not found in 'only_dep': seed off
    // (4) varname found in 'only_dep':     seed default (potentially only for selected array components)
    return dbad_json_get_seedspec("only_dep", varname, vtyp);
}

/***************************************************/
/*        end JSON stuff
/**************************************************/




static double dbad_nextRandom() {
  dbad_currentSeed += dbad_seed ;
  if (dbad_currentSeed>=1.0) dbad_currentSeed-=1.0 ;
  /* Return a value in range [1.0 2.0[ */
  return dbad_currentSeed+1.0 ;
}

static double dbad_nextIndep() {
    return dbad_nextRandom();
}

static double dbad_nextDep() {
    return dbad_nextRandom();
}

static void seed_spec_setReal8(const char* varname, seed_spec_t * seedspec, double* value) {
    if( seedspec==NULL ) {// default
	*value = dbad_nextRandom();
    }
    else if (seedspec->ndir==SEED_DEFAULT) {
        *value = dbad_nextRandom();
    }
    else if ((seedspec->ndir)==SEED_OFF) {
	*value = 0.;
    }
    else {
        fprintf(stderr, "ERROR:: ***%s*** inconsistent seed specification for variable -->%s<--\n",
                __func__, varname);
        exit(EXIT_FAILURE);
    }
}

static void seed_spec_setReal4(const char* varname, seed_spec_t * seedspec, float* value) {
    if( seedspec==NULL ) {// default
	*value = (float)dbad_nextRandom();
    }
    else if (seedspec->ndir==SEED_DEFAULT) {
        *value = (float)dbad_nextRandom();
    }
    else if ((seedspec->ndir)==SEED_OFF) {
	*value = 0.;
    }
    else {
        fprintf(stderr, "ERROR:: ***%s*** inconsistent seed specification for variable -->%s<--\n",
                __func__, varname);
        exit(EXIT_FAILURE);
    }
}

static void seed_spec_setReal8Array(const char* varname, seed_spec_t * seedspec, double* value, int length) {
    if( seedspec==NULL ) {
        for (int i=0 ; i<length ; ++i) value[i] = dbad_nextRandom();
    }
    else if( seedspec->ndir==SEED_DEFAULT ) {
        for (int i=0 ; i<length ; ++i) value[i] = dbad_nextRandom();
    }
    else if( seedspec->ndir==SEED_OFF ) {
        for (int i=0 ; i<length ; ++i) value[i] = 0.;
    }
    else if( seedspec->ndir>0 ){
        for (int i=0 ; i<length ; ++i) value[i] = 0.;
        for(int i=0; i<seedspec->ndir; ++i) {
            int j = seedspec->dir_lst[i]; // this is still 1-indexed
            if(j<1) {
                fprintf(stderr, "deteced seed index j=%d smaller than 1 (variable -->%s<--)\n",
                        j, varname);
            }
            else if( j-1>=length ) {
                fprintf(stderr, "deteced seed index j=%d larger than length=%d (variable -->%s<--)\n",
                        j, length, varname);
            }
            else {
                value[j-1] = dbad_nextRandom();
            }
        }
    }
    else {
        fprintf(stderr, "ERROR:: ***%s*** inconsistent seed specification for variable -->%s<--\n",
                __func__, varname);
        exit(EXIT_FAILURE);
    }
}

static void seed_spec_setReal4Array(const char* varname, seed_spec_t * seedspec, float* value, int length) {
    if( seedspec==NULL ) {
        for (int i=0 ; i<length ; ++i) value[i] = (float)dbad_nextRandom();
    }
    else if( seedspec->ndir==SEED_DEFAULT ) {
        for (int i=0 ; i<length ; ++i) value[i] = (float)dbad_nextRandom();
    }
    else if( seedspec->ndir==SEED_OFF ) {
        for (int i=0 ; i<length ; ++i) value[i] = 0.;
    }
    else if( seedspec->ndir>0 ){
        for (int i=0 ; i<length ; ++i) value[i] = 0.;
        for(int i=0; i<seedspec->ndir; ++i) {
            int j = seedspec->dir_lst[i]; // this is still 1-indexed
            if(j<1) {
                fprintf(stderr, "deteced seed index j=%d smaller than 1 (variable -->%s<--)\n",
                        j, varname);
            }
            else if( j-1>=length ) {
                fprintf(stderr, "deteced seed index j=%d larger than length=%d (variable -->%s<--)\n",
                        j, length, varname);
            }
            else {
                value[j-1] = dbad_nextRandom();
            }
        }
    }
    else {
        fprintf(stderr, "ERROR:: ***%s*** inconsistent seed specification for variable -->%s<--\n",
                __func__, varname);
        exit(EXIT_FAILURE);
    }
}

static void seed_spec_setComplex16(const char* varname, seed_spec_t * seedspec, double complex *value) {

    double rdot, idot;
    seed_spec_setReal8(varname, seedspec, &rdot);
    seed_spec_setReal8(varname, seedspec, &idot);
    *value = rdot + I*idot ;
}

static void seed_spec_setComplex16Array(const char* varname, seed_spec_t * seedspec, double complex *value, int length) {

    double rdot[length], idot[length];
    seed_spec_setReal8Array(varname, seedspec, rdot, length);
    seed_spec_setReal8Array(varname, seedspec, idot, length);
    for(int i=0; i<length; ++i) {
        value[i] = rdot[i] + I*idot[i] ;
    }
}

static void seed_spec_setComplex8(const char* varname, seed_spec_t * seedspec, ccmplx *value) {

    float rdot, idot;
    seed_spec_setReal4(varname, seedspec, &rdot);
    seed_spec_setReal4(varname, seedspec, &idot);
    value->r = rdot;
    value->i = idot ;
}

static void seed_spec_setComplex8Array(const char* varname, seed_spec_t * seedspec, ccmplx *value, int length) {

    float rdot[length], idot[length];
    seed_spec_setReal4Array(varname, seedspec, rdot, length);
    seed_spec_setReal4Array(varname, seedspec, idot, length);
    for(int i=0; i<length; ++i) {
        value[i].r = rdot[i];
        value[i].i = idot[i];
    }
}

void dbad_set_phase(dbad_mode_t mode) {
    char* phase = getenv("DBAD_PHASE") ;
    if(mode==ADJOINT) {
        if (phase==NULL) {
            dbad_phase = 0 ;
        } else if (strcmp(phase,"99")==0) {
            dbad_phase = 99 ;
        } else {
            dbad_phase = 0 ;
        }
    }
    else {//TANGENT
        if (phase==NULL) {
            printf("Please set DBAD_PHASE environment variable to 1 (perturbed) or 2 (tangent)\n") ;
            exit(0) ;
        } else if (strcmp(phase,"2")==0) {
            dbad_phase = 2 ;
        } else if (strcmp(phase,"1")==0) {
            dbad_phase = 1 ;
        } else if (strcmp(phase,"99")==0) {
            dbad_phase = 99 ;
        } else if (strcmp(phase,"999")==0) {
            dbad_phase = 999 ;
        } else {
            printf("DBAD_PHASE environment variable must be set to 1 or 2\n") ;
            exit(0) ;
        }
    }
}


void dbad_terminate() {
    // dispose memory
    glob_json_obj = NULL;
    json_free(&glob_json_element);
}

/*********************************************************************/
/*                                                                   */
/*        p u b l i c   i n t e r f a c e                            */
/*                                                                   */
/*********************************************************************/
void adContextTgt_init(double epsilon, double seed) {
    // set phase
    dbad_set_phase(TANGENT);

    // default settings
    dbad_mode = 1 ;
    dbad_ddeps = epsilon ;
    dbad_seed = seed ;

    // potentially override with settings from json file
    dbad_json_load();

    dbad_json_set_ddeps();
    dbad_json_set_seed();

    if (dbad_phase==2) {
        printf("Tangent code,  seed=%7.1e\n", dbad_seed) ;
        printf("=============================================\n") ;
        dbad_currentSeed = 0.0 ;
    } else if (dbad_phase==1) {
        printf("Perturbed run, seed=%7.1e, epsilon=%7.1e\n", dbad_seed, dbad_ddeps) ;
        printf("=============================================\n") ;
        dbad_currentSeed = 0.0 ;
    } else if (dbad_phase>=99) {
        printf("INTERNAL INTERFACE TESTS, seed=%7.1e, epsilon=%7.1e\n", dbad_seed, dbad_ddeps) ;
        printf("=============================================\n") ;
    } else {
        printf("DBAD_PHASE environment variable must be set to 1 or 2\n") ;
        exit(0) ;
    }
}

void adContextTgt_initReal8(char* varname, double *indep, double *indepd) {

    // look for seed specification
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    // set value
    seed_spec_setReal8(varname, seedspec, indepd);

    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase==1)
    *indep = (*indep)+dbad_ddeps*(*indepd) ;
  else if (dbad_phase>=99)
      printf("initReal8 of %s: %24.16e //%24.16e\n", varname, *indep, *indepd) ;
}

void adContextTgt_initReal8Array(char* varname, double *indep, double *indepd, int length) {
  int i ;

  // look for seed specification
  seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);

  seed_spec_setReal8Array(varname, seedspec, indepd, length);

  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i]+dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase>=99) {
    printf("initReal8Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e //%24.16e",i,indep[i],indepd[i]) ;
    printf("\n") ;
  }
}

void adContextTgt_initReal4(char* varname, float *indep, float *indepd) {
    // look for seed specification
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);

    // set value
    seed_spec_setReal4(varname, seedspec, indepd);

    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase==1)
    *indep = (*indep)+dbad_ddeps*(*indepd) ;
  else if (dbad_phase>=99)
    printf("initReal4 of %s: %24.16e //%24.16e\n", varname, *indep, *indepd) ;
}

void adContextTgt_initReal4Array(char* varname, float *indep, float *indepd, int length) {
  int i ;

  seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);

  seed_spec_setReal4Array(varname, seedspec, indepd, length);

  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i]+dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase>=99) {
    printf("initReal4Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e //%24.16e",i,indep[i],indepd[i]) ;
    printf("\n") ;
  }
}

void adContextTgt_initComplex16(char* varname, double complex *indep, double complex *indepd) {
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    seed_spec_setComplex16(varname, seedspec, indepd);
    if(seedspec!=NULL) seed_spec_dispose(seedspec);

  if (dbad_phase==1) {
    *indep = *indep + dbad_ddeps*(*indepd) ;
  } else if (dbad_phase>=99)
    printf("initComplex16 of %s: %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, creal(*indep), cimag(*indep), creal(*indepd), cimag(*indepd)) ;
}

void adContextTgt_initComplex16Array(char* varname, double complex *indep, double complex *indepd, int length) {
  int i ;
  seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);
  seed_spec_setComplex16Array(varname, seedspec, indepd, length);
  if(seedspec!=NULL) seed_spec_dispose(seedspec);

  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i] + dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase>=99) {
    printf("initComplex16Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e //%24.16e+i%24.16e",
             i,creal(indep[i]),cimag(indep[i]),creal(indepd[i]),cimag(indepd[i])) ;
    printf("\n") ;
  }
}

void adContextTgt_initComplex8(char* varname, ccmplx *indep, ccmplx *indepd) {
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    seed_spec_setComplex8(varname, seedspec, indepd);
    if(seedspec!=NULL) seed_spec_dispose(seedspec);
    
  if (dbad_phase==1) {
    indep->r = indep->r + dbad_ddeps*indepd->r ;
    indep->i = indep->i + dbad_ddeps*indepd->i ;
  } else if (dbad_phase>=99)
    printf("initComplex8 of %s: %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, indep->r, indep->i, indepd->r, indepd->i) ;
}

void adContextTgt_initComplex8Array(char* varname, ccmplx *indep, ccmplx *indepd, int length) {
  int i ;
  seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);
  seed_spec_setComplex8Array(varname, seedspec, indepd, length);
  if(seedspec!=NULL) seed_spec_dispose(seedspec);

  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i].r = indep[i].r+dbad_ddeps*indepd[i].r ;
      indep[i].i = indep[i].i+dbad_ddeps*indepd[i].i ;
    }
  } else if (dbad_phase>=99) {
    printf("initComplex8Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e //%24.16e+i%24.16e",
             i,indep[i].r,indep[i].i,indepd[i].r,indepd[i].i) ;
    printf("\n") ;
  }
}

void adContextTgt_startConclude() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_val = 0.0 ;
  dbad_condensed_tgt = 0.0 ;
}

void adContextTgt_concludeReal8(char* varname, double dep, double depd) {
    double depb;

    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);

    // set value
    seed_spec_setReal8(varname, seedspec, &depb);

    if( seedspec!=NULL ) seed_spec_dispose(seedspec);


  dbad_condensed_val += depb*(dep) ;
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depb*(depd) ;
  else if (dbad_phase>=99)
    printf("concludeReal8 of %s [%24.16e *] %24.16e //%24.16e\n", varname, depb, dep, depd) ;
}

void adContextTgt_concludeReal8Array(char* varname, double *dep, double *depd, int length) {
  int i ;
  double depb[length];

  if (dbad_phase>=99) printf("%s of %s, length=%i:\n", __func__, varname, length) ;

  seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);

  seed_spec_setReal8Array(varname, seedspec, depb, length);

  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  for (i=0 ; i<length ; ++i) {
      dbad_condensed_val += depb[i]*dep[i] ;
      if (dbad_phase==2 || dbad_phase==1) {
	  dbad_condensed_tgt += depb[i]*depd[i] ;
      } else if (dbad_phase>=99) {
	  printf("    %i:[%24.16e *] %24.16e //%24.16e",i,depb[i],dep[i],depd[i]) ;
      }
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextTgt_concludeReal4(char* varname, float dep, float depd) {
    float depb;

    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);

    seed_spec_setReal4(varname, seedspec, &depb);

    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  dbad_condensed_val += depb*(dep) ;
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depb*(depd) ;
  else if (dbad_phase>=99)
    printf("concludeReal4 of %s [%24.16e *] %24.16e //%24.16e\n", varname, depb, dep, depd) ;
}

void adContextTgt_concludeReal4Array(char* varname, float *dep, float *depd, int length) {
  int i ;
  float depb[length] ;
  if (dbad_phase>=99) printf("concludeReal8Array of %s, length=%i:\n", varname, length) ;

  seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);

  seed_spec_setReal4Array(varname, seedspec, depb, length);
  
  if (dbad_phase>=99) printf("concludeReal4Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    dbad_condensed_val += depb[i]*dep[i] ;
    if (dbad_phase==2 || dbad_phase==1) {
       dbad_condensed_tgt += depb[i]*depd[i] ;
    } else if (dbad_phase>=99) {
      printf("    %i:[%24.16e *] %24.16e //%24.16e",i,depb[i],dep[i],depd[i]) ;
    }
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextTgt_concludeComplex16(char* varname, double complex dep, double complex depd) {
    double depbr, depbi;
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);
    seed_spec_setReal8(varname, seedspec, &depbr);
    seed_spec_setReal8(varname, seedspec, &depbi);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);
    
  dbad_condensed_val += depbr*creal(dep) + depbi*cimag(dep);
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depbr*creal(depd) + depbi*cimag(depd) ;
  else if (dbad_phase>=99)
    printf("concludeComplex16 of %s [%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, depbr, depbi, creal(dep), cimag(dep), creal(depd), cimag(depd)) ;
}

void adContextTgt_concludeComplex16Array(char* varname, double complex *dep, double complex *depd, int length) {
  int i ;
  double depbr[length], depbi[length];
  seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);
  seed_spec_setReal8Array(varname, seedspec, depbr, length);
  seed_spec_setReal8Array(varname, seedspec, depbi, length);
  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99) printf("concludeComplex16Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    dbad_condensed_val += depbr[i]*creal(dep[i]) + depbi[i]*cimag(dep[i]);
    if (dbad_phase==2 || dbad_phase==1) {
      dbad_condensed_tgt += depbr[i]*creal(depd[i]) + depbi[i]*cimag(depd[i]) ;
    } else if (dbad_phase>=99) {
      printf("    %i:[%24.16e;%24.16e *] %24.16e //%24.16e",
             i, depbr[i], depbi[i], creal(dep[i]), cimag(dep[i]), creal(depd[i]), cimag(depd[i])) ;
    }
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextTgt_concludeComplex8(char* varname, ccmplx *dep, ccmplx *depd) {
    float depbr, depbi;
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);
    seed_spec_setReal4(varname, seedspec, &depbr);
    seed_spec_setReal4(varname, seedspec, &depbi);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  dbad_condensed_val += depbr*(dep->r) + depbi*(dep->i) ;
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depbr*(depd->r) + depbi*(depd->i) ;
  else if (dbad_phase>=99)
    printf("concludeComplex8 of %s [%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, depbr, depbi, dep->r, dep->i, depd->r, depd->i) ;
}

void adContextTgt_concludeComplex8Array(char* varname, ccmplx *dep, ccmplx *depd, int length) {
  int i ;
  float depbr[length], depbi[length] ;
  seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);
  seed_spec_setReal4Array(varname, seedspec, depbr, length);
  seed_spec_setReal4Array(varname, seedspec, depbi, length);
  if( seedspec!=NULL ) seed_spec_dispose(seedspec);
  if (dbad_phase>=99) printf("concludeComplex8Array of %s, length=%i:\n", varname, length) ;

  for (i=0 ; i<length ; ++i) {
    dbad_condensed_val += depbr[i]*(dep[i].r) + depbi[i]*(dep[i].i) ;
    if (dbad_phase==2 || dbad_phase==1) {
      dbad_condensed_tgt += depbr[i]*(depd[i].r) + depbi[i]*(depd[i].i) ;
    } else if (dbad_phase>=99) {
      printf("    %i:[%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e",
             i, depbr[i], depbi[i], dep[i].r, dep[i].i, depd[i].r, depd[i].i) ;
    }
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextTgt_conclude() {
  if (dbad_phase==2) {
    printf("[seed:%7.1e] Condensed result : %24.16e\n", dbad_seed, dbad_condensed_val) ;
    printf("[seed:%7.1e] Condensed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  } else if (dbad_phase==1) {
    printf("[seed:%7.1e] Condensed perturbed result : %24.16e (epsilon:%7.1e)\n",
           dbad_seed, dbad_condensed_val, dbad_ddeps) ;
    printf("[seed:%7.1e] Condensed perturbed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  }

  // final cleanup
  dbad_terminate();
}

void adContextAdj_init(double seed) {
    // set phase
    dbad_set_phase(ADJOINT);

    // default settings
    dbad_mode = 0 ;
    dbad_seed = seed ;

    // potentially override with settings from json file
    dbad_json_load();

    dbad_json_set_ddeps();
    dbad_json_set_seed();

//    char* phase = getenv("DBAD_PHASE") ;
    if(dbad_phase>=99) {
        printf("INTERNAL INTERFACE TESTS, seed=%7.1e\n", seed) ;
    }
    printf("Adjoint code,  seed=%7.1e\n", seed) ;
    printf("===================================\n") ;
    dbad_currentSeed = 0.0 ;
}

void adContextAdj_initReal8(char* varname, double *dep, double *depb) {
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);
    seed_spec_setReal8(varname, seedspec, depb);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);
  if (dbad_phase>=99)
    printf("initReal8 of %s %24.16e\n", varname, *depb) ;
}

void adContextAdj_initReal8Array(char* varname, double *dep, double *depb, int length) {
  int i ;
  seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);
  seed_spec_setReal8Array(varname, seedspec, depb, length);
  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99) {
    printf("initReal8Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e", i, depb[i]) ;
    printf("\n") ;
  }
}

void adContextAdj_initReal4(char* varname, float *dep, float *depb) {
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);
    seed_spec_setReal4(varname, seedspec, depb);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99)
    printf("initReal4 of %s %24.16e\n", varname, *depb) ;
}

void adContextAdj_initReal4Array(char* varname, float *dep, float *depb, int length) {
  int i ;
  seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);
  seed_spec_setReal4Array(varname, seedspec, depb, length);
  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99) {
    printf("initReal4Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e",i, depb[i]) ;
    printf("\n") ;
  }
}

void adContextAdj_initComplex16(char* varname, double complex *dep, double complex *depb) {
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);
    seed_spec_setComplex16(varname, seedspec, depb);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99)
    printf("initComplex16 of %s %24.16e+i%24.16e\n", varname, creal(*depb), cimag(*depb)) ;
}

void adContextAdj_initComplex16Array(char* varname, double complex *dep, double complex *depb, int length) {
  int i ;
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);
    seed_spec_setComplex16Array(varname, seedspec, depb, length);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99) {
    printf("initComplex16Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e",i, creal(depb[i]), cimag(depb[i])) ;
    printf("\n") ;
  }
}

void adContextAdj_initComplex8(char* varname, ccmplx *dep, ccmplx *depb) {
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_SCALAR);
    seed_spec_setComplex8(varname, seedspec, depb);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);
    
  if (dbad_phase>=99)
    printf("initComplex8 of %s %24.16e+i%24.16e\n", varname, depb->r, depb->i) ;
}

void adContextAdj_initComplex8Array(char* varname, ccmplx *dep, ccmplx *depb, int length) {
  int i ;
    seed_spec_t * seedspec = dbad_json_get_dep(varname, VAR_ARRAY);
    seed_spec_setComplex8Array(varname, seedspec, depb, length);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99) {
    printf("initComplex8Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e", i, depb[i].r, depb[i].i) ;
    printf("\n") ;
  }
}

void adContextAdj_startConclude() {
  dbad_currentSeed = 0.0 ;
  dbad_condensed_adj = 0.0 ;
}

void adContextAdj_concludeReal8(char* varname, double dep, double depb) {
    double depd;
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    seed_spec_setReal8(varname, seedspec, &depd);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  dbad_condensed_adj += depd*depb ;
  if (dbad_phase>=99)
    printf("concludeReal8 of %s [%24.16e *]%24.16e\n", varname, depd, depb) ;
}

void adContextAdj_concludeReal8Array(char* varname, double *dep, double *depb, int length) {
  int i ;
  double depd[length] ;
  if (dbad_phase>=99) printf("concludeReal8Array of %s, length=%i:\n", varname, length) ;
  
  seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);
  seed_spec_setReal8Array(varname, seedspec, depd, length);
  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  for (i=0 ; i<length ; ++i) {
    dbad_condensed_adj += depd[i]*depb[i] ;
    if (dbad_phase>=99) printf("    %i:[%24.16e *] %24.16e",i,depd[i],depb[i]) ;
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextAdj_concludeReal4(char* varname, float dep, float depb) {
    float depd;
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    seed_spec_setReal4(varname, seedspec, &depd);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  dbad_condensed_adj += depd*depb ;
  if (dbad_phase>=99)
    printf("concludeReal4 of %s [%24.16e *]%24.16e\n", varname, depd, depb) ;
}

void adContextAdj_concludeReal4Array(char* varname, float *dep, float *depb, int length) {
  int i ;
  float depd[length] ;
  if (dbad_phase>=99) printf("concludeReal4Array of %s, length=%i:\n", varname, length) ;
  seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);
  seed_spec_setReal4Array(varname, seedspec, depd, length);
  if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  for (i=0 ; i<length ; ++i) {
    dbad_condensed_adj += depd[i]*depb[i] ;
    if (dbad_phase>=99) printf("    %i:[%24.16e *] %24.16e",i,depd,depb[i]) ;
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextAdj_concludeComplex16(char* varname, double complex dep, double complex depb) {
    double depdr, depdi;
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    seed_spec_setReal8(varname, seedspec, &depdr);
    seed_spec_setReal8(varname, seedspec, &depdi);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  dbad_condensed_adj += depdr*creal(depb) + depdi*cimag(depb) ;
  if (dbad_phase>=99)
    printf("concludeComplex16 of %s [%24.16e+i%24.16e *]%24.16e+i%24.16e\n", varname, depdr, depdi, creal(depb), cimag(depb)) ;
}

void adContextAdj_concludeComplex16Array(char* varname, double complex *dep, double complex *depb, int length) {
  int i ;
    double depdr[length], depdi[length];
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);
    seed_spec_setReal8Array(varname, seedspec, depdr, length);
    seed_spec_setReal8Array(varname, seedspec, depdi, length);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);
    
  if (dbad_phase>=99) printf("concludeComplex16Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    dbad_condensed_adj += depdr[i]*creal(depb[i]) + depdi[i]*cimag(depb[i]) ;
    if (dbad_phase>=99) printf("    %i:[%24.16e+i%24.16e *] %24.16e+i%24.16e",i,depdr[i],depdi[i],creal(depb[i]),cimag(depb[i])) ;
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextAdj_concludeComplex8(char* varname, ccmplx *dep, ccmplx *depb) {
    float depdr, depdi;
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_SCALAR);
    seed_spec_setReal4(varname, seedspec, &depdr);
    seed_spec_setReal4(varname, seedspec, &depdi);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  dbad_condensed_adj += depdr*depb->r + depdi*depb->i ;
  if (dbad_phase>=99)
    printf("concludeComplex8 of %s [%24.16e+i%24.16e *]%24.16e+i%24.16e\n", varname, depdr, depdi, depb->r, depb->i) ;
}

void adContextAdj_concludeComplex8Array(char* varname, ccmplx *dep, ccmplx *depb, int length) {
  int i ;
    float depdr[length], depdi[length];
    seed_spec_t * seedspec = dbad_json_get_indep(varname, VAR_ARRAY);
    seed_spec_setReal4Array(varname, seedspec, depdr, length);
    seed_spec_setReal4Array(varname, seedspec, depdi, length);
    if( seedspec!=NULL ) seed_spec_dispose(seedspec);

  if (dbad_phase>=99) printf("concludeComplex8Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    dbad_condensed_adj += depdr[i]*depb[i].r + depdi[i]*depb[i].i ;
    if (dbad_phase>=99) printf("    %i:[%24.16e+i%24.16e *] %24.16e+i%24.16e",i,depdr[i],depdi[i],depb[i].r,depb[i].i) ;
  }
  if (dbad_phase>=99) printf("\n") ;
}

void adContextAdj_conclude() {
  printf("[seed:%7.1e] Condensed adjoint: %24.16e\n", dbad_seed, dbad_condensed_adj) ;
  // final cleanup
  dbad_terminate();
}

//############## INTERFACE PROCEDURES CALLED FROM FORTRAN ################

void adcontexttgt_init_(double *epsilon, double *seed) {
  adContextTgt_init(*epsilon, *seed) ;
}

void adcontexttgt_initreal8_(char* varname, double *indep, double *indepd) {
  adContextTgt_initReal8(varname, indep, indepd) ;
}

void adcontexttgt_initreal8array_(char* varname, double *indep, double *indepd, int *length) {
  adContextTgt_initReal8Array(varname, indep, indepd, *length) ;
}

void adcontexttgt_initreal4_(char* varname, float *indep, float *indepd) {
  adContextTgt_initReal4(varname, indep, indepd) ;
}

void adcontexttgt_initreal4array_(char* varname, float *indep, float *indepd, int *length) {
  adContextTgt_initReal4Array(varname, indep, indepd, *length) ;
}

void adcontexttgt_initcomplex16_(char* varname, cdcmplx *indep, cdcmplx *indepd) {
  adContextTgt_initComplex16(varname, (double complex *)indep, (double complex *)indepd) ;
}

void adcontexttgt_initcomplex16array_(char* varname, cdcmplx *indep, cdcmplx *indepd, int *length) {
  adContextTgt_initComplex16Array(varname, (double complex *)indep, (double complex *)indepd, *length) ;
}

void adcontexttgt_initcomplex8_(char* varname, ccmplx *indep, ccmplx *indepd) {
  adContextTgt_initComplex8(varname, indep, indepd) ;
}

void adcontexttgt_initcomplex8array_(char* varname, ccmplx *indep, ccmplx *indepd, int *length) {
  adContextTgt_initComplex8Array(varname, indep, indepd, *length) ;
}

void adcontexttgt_startconclude_() {
  adContextTgt_startConclude() ;
}

void adcontexttgt_concludereal8_(char* varname, double *dep, double *depd) {
  if (dbad_phase>=99)
      printf("concludereal8_ of %s: \n", varname);
  adContextTgt_concludeReal8(varname, *dep, *depd) ;
}

void adcontexttgt_concludereal8array_(char* varname, double *dep, double *depd, int *length) {
  if (dbad_phase>=99)
      printf("concludereal8array_ of %s: \n", varname);
  adContextTgt_concludeReal8Array(varname, dep, depd, *length) ;
}

void adcontexttgt_concludereal4_(char* varname, float *dep, float *depd) {
  adContextTgt_concludeReal4(varname, *dep, *depd) ;
}

void adcontexttgt_concludereal4array_(char* varname, float *dep, float *depd, int *length) {
  adContextTgt_concludeReal4Array(varname, dep, depd, *length) ;
}

void adcontexttgt_concludecomplex16_(char* varname, cdcmplx *dep, cdcmplx *depd) {
  adContextTgt_concludeComplex16(varname, *((double complex *)dep), *((double complex *)depd)) ;
}

void adcontexttgt_concludecomplex16array_(char* varname, cdcmplx *dep, cdcmplx *depd, int *length) {
  adContextTgt_concludeComplex16Array(varname, (double complex *)dep, (double complex *)depd, *length) ;
}

void adcontexttgt_concludecomplex8_(char* varname, ccmplx *dep, ccmplx *depd) {
  if (dbad_phase>=99)
      printf("concludecomplex8_ of %s: \n", varname);
  adContextTgt_concludeComplex8(varname, dep, depd) ;
}

void adcontexttgt_concludecomplex8array_(char* varname, ccmplx *dep, ccmplx *depd, int *length) {
  if (dbad_phase>=99)
      printf("concludecomplex8array_ of %s: \n", varname);
  adContextTgt_concludeComplex8Array(varname, dep, depd, *length) ;
}

void adcontexttgt_conclude_() {
  adContextTgt_conclude() ;
}

void adcontextadj_init_(double *seed) {
  adContextAdj_init(*seed) ;
}

void adcontextadj_initreal8_(char* varname, double *dep, double *depb) {
  if (dbad_phase>=99)
    printf("initreal8_ of %s \n", varname) ;
  adContextAdj_initReal8(varname, dep, depb) ;
}

void adcontextadj_initreal8array_(char* varname, double *dep, double *depb, int *length) {
  if (dbad_phase>=99)
    printf("initreal8array_ of %s \n", varname) ;
  adContextAdj_initReal8Array(varname, dep, depb, *length) ;
}

void adcontextadj_initreal4_(char* varname, float *dep, float *depb) {
  adContextAdj_initReal4(varname, dep, depb) ;
}

void adcontextadj_initreal4array_(char* varname, float *dep, float *depb, int *length) {
  adContextAdj_initReal4Array(varname, dep, depb, *length) ;
}

void adcontextadj_initcomplex16_(char* varname, cdcmplx *dep, cdcmplx *depb) {
  adContextAdj_initComplex16(varname, (double complex *)dep, (double complex *)depb) ;
}

void adcontextadj_initcomplex16array_(char* varname, cdcmplx *dep, cdcmplx *depb, int *length) {
  adContextAdj_initComplex16Array(varname, (double complex *)dep, (double complex *)depb, *length) ;
}

void adcontextadj_initcomplex8_(char* varname, ccmplx *dep, ccmplx *depb) {
  adContextAdj_initComplex8(varname, dep, depb) ;
}

void adcontextadj_initcomplex8array_(char* varname, ccmplx *dep, ccmplx *depb, int *length) {
  adContextAdj_initComplex8Array(varname, dep, depb, *length) ;
}

void adcontextadj_startconclude_() {
  adContextAdj_startConclude() ;
}

void adcontextadj_concludereal8_(char* varname, double *dep, double *depb) {
  if (dbad_phase>=99)
    printf("concludereal8_ of %s \n", varname) ;
  adContextAdj_concludeReal8(varname, *dep, *depb) ;
}

void adcontextadj_concludereal8array_(char* varname, double *dep, double *depb, int *length) {
  if (dbad_phase>=99)
    printf("concludereal8array_ of %s \n", varname) ;
  adContextAdj_concludeReal8Array(varname, dep, depb, *length) ;
}

void adcontextadj_concludereal4_(char* varname, float *dep, float *depb) {
  if (dbad_phase>=99)
    printf("concludereal4_ of %s \n", varname) ;
  adContextAdj_concludeReal4(varname, *dep, *depb) ;
}

void adcontextadj_concludereal4array_(char* varname, float *dep, float *depb, int *length) {
  if (dbad_phase>=99)
    printf("concludereal4array_ of %s \n", varname) ;
  adContextAdj_concludeReal4Array(varname, dep, depb, *length) ;
}

void adcontextadj_concludecomplex16_(char* varname, cdcmplx *dep, cdcmplx *depb) {
  adContextAdj_concludeComplex16(varname, *((double complex *)dep), *((double complex *)depb)) ;
}

void adcontextadj_concludecomplex16array_(char* varname, cdcmplx *dep, cdcmplx *depb, int *length) {
  adContextAdj_concludeComplex16Array(varname, (double complex *)dep, (double complex *)depb, *length) ;
}

void adcontextadj_concludecomplex8_(char* varname, ccmplx *dep, ccmplx *depb) {
  adContextAdj_concludeComplex8(varname, dep, depb) ;
}

void adcontextadj_concludecomplex8array_(char* varname, ccmplx *dep, ccmplx *depb, int *length) {
  adContextAdj_concludeComplex8Array(varname, dep, depb, *length) ;
}

void adcontextadj_conclude_() {
  adContextAdj_conclude() ;
}
