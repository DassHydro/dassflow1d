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
#   Model and Compilation Options
#======================================================================================================================#
COMPILO = 0

OPTIM = 0

ADJOINT = 0

SOLVER = 1

MINMETHOD = 1

VALID = 0

REF_HEIGHT = 0

# CUSTOM_DEFINES = -DREGUL_D1_ALPHA -DREGUL_BATHY
