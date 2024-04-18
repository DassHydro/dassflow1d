#======================================================================================================================!
#
#                    DassFlow1D Version 2.1
#
#======================================================================================================================!
#
#  Copyright University of Toulouse-INSA & CNRS (France)
#
#  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
#  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
#  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
#  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
#
#  DassFlow software includes few mostly independent "modules" with common architectures and structures:
#    - 1D Shallow Module (1D Shallow Water Model, Finite Volume/Finite Difference Methods), i.e. the present code.
#    - 2D Shallow Module (2D Shallow Water Model, Finite Volume Method)
#    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
#  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
#
#  Many people have contributed to the DassFlow development from the initial version to the latest ones.
#  Current main developers or scientific contributers are:
#               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
#               K. Larnier (C.S Communication and Systems & INSA Toulouse).
#               J. Monnier (INSA Toulouse & Mathematics Institute of Toulouse IMT).
#               J.-P. Vila (INSA Toulouse & Mathematics Institute of Toulouse IMT).
#               P.-A. Garambois (INSA Strasbourg & ICUBE).
#               L. Pujol (CNES & INSA Strasbourg & ICUBE).
#  and former other developers (P. Brisset, R. Madec, M. Honnorat and J. Marin).
#
#  Scientific Contact : jerome.monnier@insa-toulouse.fr
#  Technical  Contact : kevin.larnier@c-s.fr
#                       frederic.couderc@math.univ-toulouse.fr
#
#  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
#  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
#  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
#  license, users are provided only with a limited warranty and the software's author, the holder of the economic
#  rights, and the successive licensors have only limited liability.
#
#  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
#  developing or reproducing the software by the user in light of its specific status of free software, that may
#  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
#  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
#  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
#  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
#  accept its terms.
#
#======================================================================================================================!

#======================================================================================================================#
#   User Options
#======================================================================================================================#
include Makefile.inc

#======================================================================================================================#
#   Compilator Flags
#======================================================================================================================#
ifeq ($(COMPILO),1)

    FC = ifort
    CC = icc

  ifeq ($(OPTIM),0)

    export CFLAGS = -fpp -module ./obj -g -traceback -check all -check noarg_temp_created -fpe0 -nofor_main -diag-disable 8291
    export FFLAGS = -nofor_main

  else

    export CFLAGS = -fpp -module ./obj -fast -nofor_main -diag-disable 8291
    export FFLAGS = -nofor_main

  endif

else

    FC = gfortran
    CC = gcc

  ifeq ($(OPTIM),0)

    export CFLAGS = -cpp -g -finit-real=nan -fbounds-check -ffree-line-length-0 -I./build/obj -J./build/obj -ffpe-trap=invalid,zero,overflow,underflow -fPIC -fallow-argument-mismatch
    export CFLAGS_USER_DATA = -g -fbounds-check -ffree-line-length-0 -I./build/obj -ffpe-trap=invalid,zero,overflow,underflow -fPIC -fallow-argument-mismatch
    export FFLAGS = -fPIC -fallow-argument-mismatch

  else

    export CFLAGS = -cpp -O3 -finit-real=nan -fno-range-check -ffree-line-length-0 -I./build/obj -J./build/obj -fPIC -fallow-argument-mismatch
    export CFLAGS_USER_DATA = -O3 -fno-range-check -ffree-line-length-0 -I./build/obj -fPIC -fallow-argument-mismatch
    export FFLAGS = -fPIC -fallow-argument-mismatch

  endif

endif

export CLIBS =
export FLIBS =

#======================================================================================================================#
#   Scotch
#======================================================================================================================#

SCOTCH = $(PWD)/libs/scotch_5.1.12_esmumps


#======================================================================================================================#
#   MUMPS
#======================================================================================================================#

MUMPS = $(PWD)/libs/MUMPS_4.10.0

ifeq ($(SOLVER),1)

  export CLIBS  +=  -I$(SCOTCH)/include

  export FLIBS  +=    $(SCOTCH)/lib/libscotch.a
  export FLIBS  +=    $(SCOTCH)/lib/libscotcherr.a

  export CLIBS  +=  -I$(MUMPS)/include -I$(MUMPS)/libseq

  export FLIBS  +=    $(MUMPS)/lib/libdmumps.a
  export FLIBS  +=    $(MUMPS)/lib/libmumps_common.a
  export FLIBS  +=  -L$(MUMPS)/PORD/lib
  export FLIBS  +=  -L$(MUMPS)/libseq
  export FLIBS  +=  -lpord -lmpiseq -lblas -lpthread -L$(PWD)/libs/scotch_5.1.12_esmumps/lib -lesmumps -lscotch -lscotcherr

endif

#======================================================================================================================#
#   AGMG
#======================================================================================================================#

AGMG = $(PWD)/libs/AGMG_3.2.0-aca/SRC

ifeq ($(SOLVER),2)

  export FLIBS  +=  $(AGMG)/dagmg.o
  export FLIBS  +=  $(AGMG)/dagmg_mumps.o
  export FLIBS  +=  -L/usr/lib
  export FLIBS  +=  -llapack -lblas

endif

#======================================================================================================================#
#   BOOST-PYTHON
#======================================================================================================================#

export CLIBS_PY  =  -I$(BOOST_PYTHON_INC) -I$(PYTHON_INC)
export CFLAGS_PY =  -O3 -fPIC
export FLIBS_PY  = $(FLIBS) $(BOOST_PYTHON_LIB) $(PYTHON_LIB)

#======================================================================================================================#
#	CPP options according the options choices above
#======================================================================================================================#

CPP_MODEL     = -DUSE_SW_MONO -DSTRICKLER_EINSTEIN -DAVERAGE_CONVEYANCE -DDEBORD_FORMULA
CPP_ADJ       =
CPP_SOL       =
CPP_VAL       =
CPP_COMP      =
CPP_REF_HEIGHT    =

ifeq ($(COMPILO),1)
  CPP_COMP = -DUSE_INTEL
endif

ifeq ($(ADJOINT),1)
  CPP_ADJ = -DUSE_ADJ
  ifeq ($(shell which tapenade),)
     TAPENADE_VERSION=
     TAPENADE_DIR=
  else
    TAPENADE_VERSION=$(shell tapenade --version | head -1 | cut -d' ' -f2)
    TAPENADE_VERSION_MAJOR=$(shell tapenade --version | head -1 | cut -d' ' -f2 | cut -d'.' -f1)
    TAPENADE_VERSION_MINOR=$(shell tapenade --version | head -1 | cut -d' ' -f2 | cut -d'.' -f2)
    TAPENADE_VERSION_PADDED=$(TAPENADE_VERSION_MAJOR).$(shell printf '%02d' $(TAPENADE_VERSION_MINOR))
    TAPENADE_DIR=$(shell tapenade --version | head -2 | tail -1 | cut -d'=' -f2)
    ifeq ($(wildcard $(TAPENADE_DIR)/ADFirstAidKit),)
        TAPENADE_DIR=$(TAPENADE_HOME)
    endif
#     TAP_INPUTVARS=dof0%q dof0%s dof0%h bc%hyd%t bc%hyd%q bc%hyd_FS%a0 bc%hyd_FS%a bc%hyd_FS%b bc%rat%h bc%rat%q K_params%alpha K_params%beta alpha_ratcurve beta_ratcurve bathy_points bathy_cell qin_chg alpha_K_chg beta_K_chg bathy_cell_chg bathy_points_chg bc%hyd_lat%t bc%hyd_lat%q
#     TAP_OUTPUTVARS=dof%q dof%s dof%h innovation%diff cost

    TAP_INPUTVARS=ctrl%x
    TAP_OUTPUTVARS=mdl%dof%q mdl%dof%h cost

    ifeq ($(shell expr $(TAPENADE_VERSION_PADDED) \< 3.13), 1)
      TAPENADE_DIFFARGS=-d -fixinterface -difffuncname _diff -diffvarname _diff -head run_model \
                      -vars "$(TAP_INPUTVARS)" -outvars "$(TAP_OUTPUTVARS)" -O ./build/tap
      TAPENADE_BACKARGS=-b -fixinterface -difffuncname _back -diffvarname _back -head run_model -html \
                      -vars "$(TAP_INPUTVARS)" -outvars "$(TAP_OUTPUTVARS)" -O ./build/tap
    else
      TAPENADE_DIFFARGS=-d -fixinterface -tgtfuncname _diff -tgtvarname _diff -tgtmodulename _diff -copyname _cd \
                      -head "calc_cost($(TAP_INPUTVARS))\($(TAP_OUTPUTVARS))" -O ./build/tap
      TAPENADE_BACKARGS=-b -fixinterface -adjfuncname _back -adjvarname _back -adjmodulename _back -copyname _cb -html \
                      -head "calc_cost($(TAP_INPUTVARS))\($(TAP_OUTPUTVARS))" -O ./build/tap
    endif
  endif
endif

ifeq ($(SOLVER),1)
  CPP_SOL = -DUSE_MUMPS
endif

ifeq ($(SOLVER),2)
  CPP_SOL = -DUSE_AGMG
endif

ifeq ($(MINMETHOD),1)
  ifeq ($(ADJOINT),1)
    CPP_MIN = -DUSE_M1QN3
    OBJS_MINMETHOD=obj/m1qn3.o\
                   obj/ddot.o
  endif
endif

ifeq ($(MINMETHOD),2)
  CPP_MIN = -DUSE_N2QN1
  OBJS_MINMETHOD=obj/n2qn1.o
endif

ifeq ($(MINMETHOD),3)
  CPP_MIN = -DUSE_LBFGSB3
  OBJS_MINMETHOD=obj/lbfgsb.o\
                 obj/lbfgsb_timer.o\
                 obj/lbfgsb_blas.o\
                 obj/lbfgsb_linpack.o
endif

ifeq ($(VALID),1)
  CPP_VAL = -DUSE_VALID
endif

ifeq ($(REF_HEIGHT),1)
  CPP_REF_HEIGHT = -DUSE_REF_HEIGHT
endif

# CPP_MISC=-DHIDE_MACHINE_LIMITS -DFLOODPLAIN_MODEL
CPP_MISC=-DHIDE_MACHINE_LIMITS $(CUSTOM_DEFINES)
CPP_FLAGS      =   $(CPP_COMP) $(CPP_MODEL) $(CPP_ADJ) $(CPP_SOL) $(CPP_MIN) $(CPP_VAL) $(CPP_REF_HEIGHT) $(CPP_MISC)
CPP_FLAGS_ADJ  =  $(CPP_COMP) $(CPP_MODEL) -DCPP_ADJ $(CPP_MISC)

#======================================================================================================================#
#
#======================================================================================================================#

MODELD = src/sw_mono

ADJD =

ifeq ($(ADJOINT),1)
  ADJD = src/adjoint
endif

DIRS      =  src src/common src/base $(MODELD) $(ADJD)
DIRS_INC  =  src/base/include $(MODELD)/include

POST_SRC  =  $(foreach DIR, $(DIRS)    , $(patsubst $(DIR)/%.f90,build/cpp/%.f90,$(wildcard $(DIR)/*.f90)))
POST_INC  =  $(foreach DIR, $(DIRS_INC), $(patsubst $(DIR)/%.f90,build/cpp/%.inc,$(wildcard $(DIR)/*.f90)))

VPATH  =  src:src/common:src/core:src/base:$(MODELD):$(ADJD):src/base/include:$(MODELD)/include

ifeq ($(ADJOINT),1)
  OBJS_TAP  =  obj/m_tap_vars.o\
               obj/m_adjoint.o\
               obj/m_numeric_back.o\
               obj/m_linear_solver_back.o\
               $(patsubst src/adjoint/%.f90,obj/%.o,$(wildcard src/adjoint/*.f90))\
               obj/adStack.o\
               $(OBJS_MINMETHOD) \
               $(patsubst tap/%.f90,obj/%.o,$(sort $(wildcard tap/m_*.f90)))\
               $(patsubst tap/%.f90,obj/%.o,$(wildcard tap/*.f90))
else
  OBJS_TAP  =
endif

OBJS_MOD  =  build/obj/m_common.o\
             build/obj/m_linear_algebra.o\
             build/obj/m_mesh.o\
             build/obj/m_time_screen.o\
             build/obj/m_numeric.o \
             build/obj/m_opt.o\
             build/obj/m_setup.o\
             $(patsubst $(MODELD)/%.f90,build/obj/%.o,$(wildcard $(MODELD)/m_*.f90))\
             build/obj/m_linear_solver.o
OBJS_SRC  =  build/obj/m_common.o\
             build/obj/m_linear_algebra.o\
             build/obj/m_mesh.o\
             build/obj/m_time_screen.o\
             build/obj/m_numeric.o\
             $(patsubst $(MODELD)/%.f90,build/obj/%.o,$(wildcard $(MODELD)/m_*.f90))\
             build/obj/m_linear_solver.o\
             $(patsubst src/common/%.f90,build/obj/%.o,$(wildcard src/common/*.f90))\
             $(patsubst src/base/%.f90,build/obj/%.o,$(wildcard src/base/*.f90))\
             $(patsubst $(MODELD)/%.f90,build/obj/%.o,$(wildcard $(MODELD)/*.f90))\
             build/obj/m_linear_solver.o\
             $(OBJS_TAP)\
             build/obj/main.o
             
OBJS_SRC_UNITTEST  =  obj/m_common.o\
                      obj/m_linear_algebra.o\
                      obj/m_mesh.o\
                      obj/m_time_screen.o\
                      obj/m_numeric.o\
                      $(patsubst $(MODELD)/%.f90,obj/%.o,$(wildcard $(MODELD)/m_*.f90))\
                      obj/m_linear_solver.o\
                      $(patsubst src/common/%.f90,obj/%.o,$(wildcard src/common/*.f90))\
                      $(patsubst src/base/%.f90,obj/%.o,$(wildcard src/base/*.f90))\
                      $(patsubst $(MODELD)/%.f90,obj/%.o,$(wildcard $(MODELD)/*.f90))\
                      $(OBJS_TAP)
             
OBJS_PY_SRC = $(patsubst src/python/%.f90,obj/%.mo,$(wildcard src/python/*.f90))\
              $(patsubst src/python/%.cpp,obj/%.mo,$(wildcard src/python/*.cpp))
              
# DEPRECATED
# TODO remove all references to CASEDIR
CASEDIR=bin

#======================================================================================================================#
#
#======================================================================================================================#
FILES_TAP = build/tap/m_common.f90 \
            build/tap/m_linear_algebra.f90 \
            build/tap/m_linear_solver.f90 \
            build/tap/m_mesh.f90 \
            build/tap/m_numeric.f90 \
            build/tap/m_sw_mono.f90 \
            build/tap/m_obs.f90 \
            build/tap/m_control.f90 \
            build/tap/bathy_slopes.f90 \
            build/tap/apply_bathy_field.f90 \
            build/tap/apply_strickler_fields.f90 \
            build/tap/geometry.f90 \
            build/tap/apply_control.f90 \
            build/tap/calc_cost.f90 \
            build/tap/calc_estimations.f90 \
            build/tap/reference_depths.f90 \
            build/tap/standard_step.f90 \
            build/tap/preissmann_timestep.f90 \
            build/tap/implicit_diffusive_wave.f90 \
            build/tap/steady_states_loop.f90 \
            build/tap/time_loop.f90

#======================================================================================================================#
#
#======================================================================================================================#
all: build/api/_dassflow1d.so
rebuild_all: cleanall generate_adjoint all
swot_discharge_algorithm_solver: ../swot_discharge_algorithm_dassflow1d/bin/dassflow1d_solver
doc_dev:
	@echo "================================================================================"
	@echo "  Generating developers documentation                                          *"
	@echo "================================================================================"
	@cd src ; doxygen doxygen_settings.txt ; cd ..
	@echo "================================================================================"
	@echo "  Developers documentation generated in doc/dev                                *"
	@echo "================================================================================"


#======================================================================================================================#
#  Preprocessor
#======================================================================================================================#

preproc: $(POST_SRC) $(POST_INC)

mes:
	@echo "================================================================================"
	@echo "  DassFlow compilation and link done                                           *"
	@echo "================================================================================"

#======================================================================================================================#
#
#======================================================================================================================#

# $(python_module): $(OBJS_SRC) $(OBJS_PY_SRC)
# 	@echo "================================================================================"
# 	$(FC) $(FFLAGS) -shared -o ./lib/$(python_module) ./obj/*.mo $(FLIBS_PY) -lstdc++

#======================================================================================================================#
#  Tapenade
#======================================================================================================================#

tapfiles: copy tgt adj post

tap_files: copy tgt adj post

#generate_adjoint: copy tgt adj post
ifeq ($(ADJOINT),1)
generate_adjoint: copy tgt adj post
else
generate_adjoint:
	$(error ADJOINT must be set to 1");
endif

tap_files_no_post: copy tgt adj

copy: $(FILES_TAP)

tgt: $(FILES_TAP)
	@echo "================================================================================"
# 	tapenade -d -fixinterface \
# 	-difffuncname _diff -diffvarname _diff -head run_model \
#    -vars "dof0%q dof0%s dof0%h bc%hyd%t bc%hyd%q bc%hyd_FS%a0 bc%hyd_FS%a bc%hyd_FS%b bc%rat%h bc%rat%q K_params%alpha K_params%beta alpha_ratcurve beta_ratcurve bathy_points bathy_cell qin_chg alpha_K_chg beta_K_chg bathy_points_chg bathy_cell_chg bc%hyd_lat%t bc%hyd_lat%q" \
# 	-outvars "dof%q dof%s dof%h dof%qlat innovation%diff cost" -O ./tap  $(FILES_TAP)
	tapenade $(TAPENADE_DIFFARGS) $(FILES_TAP)
adj: $(FILES_TAP)
	@echo "================================================================================"
# 	tapenade -b -fixinterface \
# 	-difffuncname _back -diffvarname _back -head run_model -html \
#    -vars "dof0%q dof0%s dof0%h bc%hyd%t bc%hyd%q bc%hyd_FS%a0 bc%hyd_FS%a bc%hyd_FS%b bc%rat%h bc%rat%q K_params%alpha K_params%beta alpha_ratcurve beta_ratcurve bathy_points bathy_cell qin_chg alpha_K_chg beta_K_chg bathy_points_chg bathy_cell_chg bc%hyd_lat%t bc%hyd_lat%q" \
# 	-outvars "dof%q dof%s dof%h dof%qlat innovation%diff cost" -O ./tap  $(FILES_TAP)
	tapenade $(TAPENADE_BACKARGS) $(FILES_TAP)

# post:
# 	@echo "====================================== ========================================="
# 	cp ./src/adjoint/finish_to_gen_adjoint.pl ./tap/
# 	cd ./tap ; perl finish_to_gen_adjoint.pl
# 	@echo "================================================================================"

post:
	@echo "====================================== ========================================="
	cp ./src/adjoint/finish_to_gen_adjoint.pl ./build/tap/
	cd ./build/tap ; perl finish_to_gen_adjoint.pl
	@echo "================================================================================"

#======================================================================================================================#
#  Unit tests
#======================================================================================================================#

unit_tests: ./src/unit_tests/m_unittest.o $(OBJS_SRC_UNITTEST) $(EXE_UNITTESTS)

./src/unit_tests/m_unittest.o: ./src/unit_tests/m_unittest.f90
	$(FC) -cpp $(CFLAGS) $(CLIBS) -o $@ -c $<

./src/unit_tests/%: ./src/unit_tests/%.f90 $(OBJS_SRC_UNITTEST) ./src/unit_tests/m_unittest.o
	$(FC) -cpp $(CFLAGS) $(CLIBS) -o $@ $^ $(FLIBS)


#======================================================================================================================#
#  Dependencies
#======================================================================================================================#

./build/cpp/%.f90: %.f90
	@echo "================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./build/cpp/%.inc: %.f90
	@echo "================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./build/tap/%.f90: %.f90
	@echo "================================================================================"
	sed -e '/<NOADJ/,/>NOADJ/d;/NOADJ/d;/^ *!.*/d' $< | cpp -C -P $(CPP_FLAGS_ADJ) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

# ./build/obj/%.o: ./build/cpp/%.f90 $(POST_INC)
# 	@echo "================================================================================"
# 	$(FC) $(CFLAGS) $(CLIBS) -c $< -o $@
./build/obj/m_linear_algebra.o: m_linear_algebra.f90 \
                      ./build/obj/m_common.o
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@

./build/obj/m_linear_solver.o: m_linear_solver.f90 \
                      ./build/obj/m_common.o\
                      ./build/obj/m_linear_algebra.o\
                      ./build/obj/m_mesh.o
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@
	
./build/obj/m_mesh.o: m_mesh.f90 \
                      ./build/obj/m_common.o \
                      ./build/obj/m_linear_algebra.o
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@
	
./build/obj/m_setup.o: m_setup.f90 \
                       ./build/obj/m_common.o \
                       ./build/obj/m_mesh.o \
                       ./build/obj/m_opt.o
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@
	
./build/obj/m_%.o: m_%.f90
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@
	
./build/obj/%.o: %.f90
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@
	
./build/obj/%.o: wrappers/%.f90
	@echo "================================================================================"
	$(FC) $(CPP_FLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@

./build/obj/%.o: ./build/tap/%.f90
	@echo "============================================================================================================"
	$(FC) $(CFLAGS) -c $< -o $@

# ./build/obj/%.o: ./tap/%.f90
# 	@echo "============================================================================================================"
# 	$(FC) $(CFLAGS) -c $< -o $@

./build/obj/admm_tapenade_interface.o: $(TAPENADE_DIR)/ADFirstAidKit/admm_tapenade_interface.f90
	@echo "================================================================================"
	$(FC) $(CFLAGS) -c $< -o $@

./build/obj/adBuffer.o: $(TAPENADE_DIR)/ADFirstAidKit/adBuffer.f
	@echo "================================================================================"
	$(FC) $(FFLAGS) -c $< -o $@

./build/obj/adStack.o: $(TAPENADE_DIR)/ADFirstAidKit/adStack.c
	@echo "================================================================================"
	$(CC) $(FFLAGS) -c $< -o $@

./build/obj/%.o: %.f
	@echo "================================================================================"
# 	$(FC) $(CFLAGS) -c $< -o $@
	$(FC) $(FFLAGS) -c $< -o $@

./build/obj/%.o: %.for
	@echo "================================================================================"
	$(FC) $(FFLAGS) -c $< -o $@

./build/obj/%.o: %.c
	@echo "================================================================================"
# 	$(CC) -O2 -c $< -o $@
	$(CC) $(CFLAGS) -c $< -o $@

./build/obj/%.mo: ./src/python/%.f90
	@echo "================================================================================"
	$(FC) -cpp -fPIC $(CFLAGS) $(CLIBS) -c $< -o $@

./build/obj/%.mo: ./src/python/%.cpp
	@echo "================================================================================"
	$(CC) $(CFLAGS_PY) $(CLIBS_PY) -c $< -o $@
	
$(CASEDIR)/m_user_data.cpp.f90: $(CASEDIR)/m_user_data.f90
	@echo "================================================================================"
	sed -e '/<NOADJ/,/>NOADJ/d;/NOADJ/d;/^ *!.*/d' $< | cpp -C -P $(CPP_FLAGS_ADJ) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@
	
$(CASEDIR)/m_user_data.o: $(CASEDIR)/m_user_data.cpp.f90
	@echo "================================================================================"
	$(FC) $(CFLAGS_USER_DATA) $(CLIBS) -c $< -o $@
	
build/obj/m_user_data.o1: src/m_user_data.f90
	@echo "================================================================================"
	$(FC) $(CFLAGS) $(CLIBS) -c $< -o $@
	
build:
	mkdir -p build

build/cpp: build
	mkdir -p build/cpp

build/obj: build
	mkdir -p build/obj

build/tap: build
	mkdir -p build/tap
	
build/api:
	mkdir -p build/api

build/api/dassflow1d: build/api
	mkdir -p build/api/dassflow1d

build/api/dassflow1d/assim:
	mkdir -p build/api/dassflow1d/assim

build/api/dassflow1d/post:
	mkdir -p build/api/dassflow1d/post

build/api/dassflow1d/utils:
	mkdir -p build/api/dassflow1d/utils

build/api/dassflow1d/assim/%.py: src/wrappers/assim/%.py build/api/dassflow1d/assim
	@echo "================================================================================"
	cp $< $@

build/api/dassflow1d/post/%.py: src/wrappers/post/%.py build/api/dassflow1d/post
	@echo "================================================================================"
	cp $< $@

build/api/dassflow1d/utils/%.py: src/wrappers/utils/%.py build/api/dassflow1d/utils
	@echo "================================================================================"
	cp $< $@
# 	
# build/api/dassflow1d/utils/__init__.py: build/api/dassflow1d/utils
# 	@echo "================================================================================"
# 	touch $@

#======================================================================================================================#
#     Libraries
#======================================================================================================================#

alllibs:
ifeq ($(SOLVER),1)
	rm -rf $(SCOTCH)/src/Makefile.inc
ifeq ($(COMPILO),0)
	cd $(SCOTCH)/src ; ln -s Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc
else
	cd $(SCOTCH)/src ; ln -s Make.inc/Makefile.inc.x86-64_pc_linux2.icc Makefile.inc
endif
	cd $(PWD)/libs/scotch_5.1.12_esmumps/src ;\
	make
	cd $(PWD)/libs/MUMPS_4.10.0 ;\
	make
endif
ifeq ($(SOLVER),2)
	cd $(AGMG) ; make dseq
endif

cleanlibs:
	cd $(PWD)/libs/scotch_5.1.12_esmumps/src ;\
	make clean ;\
	cd .. ;\
	rm -rf bin lib include ;\
	cd $(PWD)/libs/MUMPS_4.10.0 ;\
	make clean

#======================================================================================================================#
#     Python wrappers
#======================================================================================================================#
# WRAPPERS_SRC = src/common/m_mesh.f90 \
#                src/sw_mono/m_obs.f90 \
#                src/sw_mono/m_control.f90 \
#                src/sw_mono/m_opt.f90 \
#                src/sw_mono/m_setup.f90 \
#                src/sw_mono/m_model.f90 \
#                src/base/geometry.f90 \
#                src/base/read_mesh.f90 \
#                src/base/standard_step.f90 \
#                build/api/apply_control.f90 \
#                build/api/preissmann_timestep.f90 \
#                build/api/time_loop.f90


#                build/api/m_opt.f90 \
#                build/api/m_setup.f90 \

# WRAPPERS_CPP = build/api/m_linear_algebra.f90 \
#                build/api/free_surface_slopes.f90 \
#                build/api/bathy_slopes.f90 \
#                build/api/curvilinear_abscissae.f90 \
#                build/api/read_mesh.f90 \
#                build/api/rectangular_mesh.f90 \
#                build/api/write_mesh.f90 \
#                build/api/read_spatial_field.f90 \
#                build/api/m_mesh.f90 \
#                build/api/m_minimization.f90 \
#                build/api/m_obs.f90 \
#                build/api/m_sw_mono.f90 \
#                build/api/m_control.f90 \
#                build/api/geometry.f90 \
#                build/api/apply_bathy_field.f90 \
#                build/api/apply_strickler_fields.f90 \
#                build/api/apply_control.f90 \
#                build/api/calc_cost.f90 \
#                build/api/generate_observations.f90 \
#                build/api/preissmann_timestep.f90 \
#                build/api/implicit_diffusive_wave.f90 \
#                build/api/standard_step.f90 \
#                build/api/steady_state.f90 \
#                build/api/steady_states_loop.f90 \
#                build/api/time_loop.f90
# #               build/api/time_loop.f90 \
# #               build/api/m1qn3.f
# #                build/api/resample_mesh.f90 \

WRAPPERS_CPP = build/api/m_linear_algebra.f90 \
               build/api/read_mesh.f90 \
               build/api/rectangular_mesh.f90 \
               build/api/write_mesh.f90 \
               build/api/m_mesh.f90 \
               build/api/m_minimization.f90 \
               build/api/m_obs.f90 \
               build/api/m_sw_mono.f90 \
               build/api/m_control.f90 \
               build/api/calc_cost.f90
               
ifeq ($(ADJOINT),1)
	WRAPPERS_CPP += build/api/calc_cost_and_gradients.f90
endif

#           ./build/obj/m_time_screen.o \
#           ./build/obj/m_opt.o \
#           ./build/obj/m_setup.o \

OBJ_PY = ./build/obj/m_common.o \
          ./build/obj/m_linear_algebra.o \
          ./build/obj/m_stdout.o \
          ./build/obj/m_mesh.o \
          ./build/obj/m_minimization.o \
          ./build/obj/m_numeric.o \
          ./build/obj/m_obs.o \
          ./build/obj/m_sw_mono.o \
          ./build/obj/m_linear_solver.o \
          ./build/obj/m_control.o \
          ./build/obj/wrappers_utils.o \
          ./build/obj/apply_bathy_field.o \
          ./build/obj/apply_strickler_fields.o \
          ./build/obj/bathy_slopes.o \
          ./build/obj/curvilinear_abscissae.o \
          ./build/obj/geometry.o \
          ./build/obj/init_implicit.o \
          ./build/obj/read_mesh.o \
          ./build/obj/rectangular_mesh.o \
          ./build/obj/write_mesh.o \
          ./build/obj/read_spatial_field.o \
          ./build/obj/resample_mesh.o\
          ./build/obj/apply_control.o \
          ./build/obj/calc_cost.o \
          ./build/obj/calc_estimations.o \
          ./build/obj/free_surface_slopes.o \
          ./build/obj/generate_observations.o \
          ./build/obj/reference_depths.o \
          ./build/obj/preissmann_timestep.o \
          ./build/obj/implicit_diffusive_wave.o \
          ./build/obj/output_sw.o \
          ./build/obj/solve_system_mumps.o \
          ./build/obj/standard_step.o \
          ./build/obj/steady_state.o \
          ./build/obj/steady_states_loop.o  \
          ./build/obj/time_loop.o \
          ./build/obj/update_internal_counters.o \
          ./build/obj/update_probes.o
ifeq ($(ADJOINT),1)
	OBJ_PY += ./build/obj/adBuffer.o \
            ./build/obj/adStack.o \
            ./build/obj/admm_tapenade_interface.o \
            ./build/obj/m_linear_algebra_back.o \
            ./build/obj/m_linear_solver_back.o \
            ./build/obj/m_mesh_back.o \
            ./build/obj/m_obs_back.o \
            ./build/obj/m_sw_mono_back.o \
            ./build/obj/m_control_back.o \
            ./build/obj/m_numeric_back.o \
            ./build/obj/calc_cost_and_gradients.o \
            ./build/obj/init_back.o \
            ./build/obj/solve_system_mumps_back.o \
            ./build/obj/m1qn3.o \
            $(patsubst ./build/tap/%.f90, ./build/obj/%.o, $(wildcard ./build/tap/*_cb.f90)) \
            $(patsubst ./build/tap/%.f90, ./build/obj/%.o, $(wildcard ./build/tap/*_back.f90))
endif
ASSIM_PY =  ./build/api/dassflow1d/assim/__init__.py \
            ./build/api/dassflow1d/assim/m1qn3.py \
            ./build/api/dassflow1d/assim/window_assim.py
POST_PY =  ./build/api/dassflow1d/post/__init__.py \
           ./build/api/dassflow1d/post/metrics.py \
           ./build/api/dassflow1d/post/results.py \
           ./build/api/dassflow1d/post/minimization.py
UTILS_PY = ./build/api/dassflow1d/utils/__init__.py \
           ./build/api/dassflow1d/utils/hecras_geoparser.py \
           ./build/api/dassflow1d/utils/SwotObs.py \
           ./build/api/dassflow1d/utils/PepsiNetCDF.py \
           ./build/api/dassflow1d/utils/mesh_from_hecras_geo.py \
           ./build/api/dassflow1d/utils/mesh_from_swot_observations.py \
           ./build/api/dassflow1d/utils/mesh_utils.py \
           ./build/api/dassflow1d/utils/pepsi_nodes_observations.py \
           ./build/api/dassflow1d/utils/sge_rivernodes_observations.py \
           ./build/api/dassflow1d/utils/stdout_utils.py

WRAPPED_OBJ = $(patsubst ./build/obj/%.o, ../obj/%.o, $(OBJ_PY))

build/api/env.sh:
	@echo 'export PYTHONPATH=$(CURDIR)/build/api:$$PYTHONPATH' > build/api/env.sh

build/api/kind_map: src/wrappers/kind_map
	@echo "============================================================================================================"
	cp $< $@

build/api/%.f90: src/common/%.f90 build/obj/%.o build/api/kind_map
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@
	
build/api/%.f90: src/adjoint/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@
	
build/api/%.f90: src/base/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@
	
build/api/%.f90: src/sw_mono/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@
	
# build/api/m1qn3.f: src/adjoint/m1qn3.f
# 	@echo "============================================================================================================"
# 	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@


copy_python_utils: $(UTILS_PY)

build/api/_dassflow1d.so: build/api/dassflow1d \
                          build/obj \
                          build/tap \
                          $(OBJ_PY) $(WRAPPERS_CPP) $(ASSIM_PY) $(POST_PY) $(UTILS_PY)
	cd build/api ; f90wrap -m dassflow1d $(patsubst build/api/%, %, ${WRAPPERS_CPP}) -k kind_map -S 128 -P -M --shorten-routine-names
	cd build/api ; f2py-f90wrap --fcompiler=gfortran --opt="-O3" --build-dir . -c -m _dassflow1d -I../obj -L. f90wrap*.f90 ../obj/*.o $(FLIBS)
	cp ./src/wrappers/finish_to_gen_wrappers.pl ./build/api/
	cd ./build/api ; perl finish_to_gen_wrappers.pl

step_wrappers: $(OBJ_PY) $(WRAPPERS_CPP) $(ASSIM_PY) $(POST_PY) $(UTILS_PY)
	cd build/api ; f2py-f90wrap --fcompiler=gfortran --opt="-O3" --build-dir . -c -m _dassflow1d -I../obj -L. f90wrap*.f90 ../obj/*.o $(FLIBS)
	cp ./src/wrappers/finish_to_gen_wrappers.pl ./build/api/
	cd ./build/api ; perl finish_to_gen_wrappers.pl

post_wrappers:
	cp ./src/wrappers/finish_to_gen_wrappers.pl ./build/api/
	cd ./build/api ; perl finish_to_gen_wrappers.pl


#======================================================================================================================#
#     Clean
#======================================================================================================================#

clean:
	@echo "================================================================================"
	@echo "  Clean compilation files                                                      *"
	@echo "================================================================================"
	@rm -rf build/obj/*
	@rm -rf build/api/*.f90

cleantap:
	@echo "================================================================================"
	@echo "*  Cleaning tap directory                                                      *"
	@echo "================================================================================"
	@rm -rf ./build/tap/*

cleanall: clean cleantap

#======================================================================================================================#
#     Misc
#======================================================================================================================#
print_params:
	@echo "COMPILO=$(COMPILO)"
	@echo "ADJOINT=$(ADJOINT)"
	@echo "TAPENADE_VERSION=$(TAPENADE_VERSION) (padded:$(TAPENADE_VERSION_PADDED))" 
	@echo "TAPENADE_DIR=$(TAPENADE_DIR)"
	@echo TAPENADE_BACKARGS=$(TAPENADE_BACKARGS)
