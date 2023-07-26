!======================================================================================================================!
!
!                    DassFlow1D Version 2.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures:
!    - 1D Shallow Module (1D Shallow Water Model, Finite Volume/Finite Difference Methods), i.e. the present code.
!    - 2D Shallow Module (2D Shallow Water Model, Finite Volume Method)
!    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current main developers or scientific contributers are:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               K. Larnier (C.S Communication and Systems & INSA Toulouse).
!               J. Monnier (INSA Toulouse & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA Toulouse & Mathematics Institute of Toulouse IMT).
!               P.-A. Garambois (INSA Strasbourg & ICUBE).
!               L. Pujol (CNES & INSA Strasbourg & ICUBE).
!  and former other developers (P. Brisset, R. Madec, M. Honnorat and J. Marin).
!
!  Scientific Contact : jerome.monnier@insa-toulouse.fr
!  Technical  Contact : kevin.larnier@c-s.fr
!                       frederic.couderc@math.univ-toulouse.fr
!
!  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
!  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
!  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
!
!  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
!  license, users are provided only with a limited warranty and the software's author, the holder of the economic
!  rights, and the successive licensors have only limited liability.
!
!  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
!  developing or reproducing the software by the user in light of its specific status of free software, that may
!  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
!  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
!  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
!  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
!
!  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
!  accept its terms.
!
!======================================================================================================================!
!> @file solve_system_mumps.f90
!> @brief This file includes the subroutine for solving the Shallow-Water Equations using Preissmann scheme with LPI 
!> and double sweep methods.


!> \brief Solve the Shallow-Water Equations using Preissmann scheme with LPI and double sweep methods.
subroutine solve_system_mumps_back(mumps_verbose, imp, imp_back)
    use m_common
    use m_sw_mono, only: ImplicitMatrix
    use m_sw_mono_back, only: ImplicitMatrix_back
    use m_numeric
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Verbose code
    integer(ip) :: mumps_verbose
    !> @param imp Implicit matrix
    type(ImplicitMatrix), intent(inout) :: imp
    !> @param imp Implicit matrix
    type(ImplicitMatrix_back), intent(inout) :: imp_back
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of column
    integer(ip) :: col
    ! Iterator
    integer(ip) :: i
    ! Index of row
    integer(ip) :: row

    imp%sys%JOB = 6
    imp%sys%ICNTL(9) = 2
!     if (mumps_verbose > 0) then
!         imp%sys%ICNTL(1) = 1
!         imp%sys%ICNTL(2) = 1
!         imp%sys%ICNTL(3) = 1
!     else
!         imp%sys%ICNTL(1) = 0
!     end if
    
    imp%sys%A(:) = imp%ANZ(:)
    imp%sys%RHS(:) = imp_back%RHS(:)
    
!     imp%sys%JOB=1
!     imp%sys%WRITE_PROBLEM="matrix1d_back"
!     call DMUMPS(imp%sys)
!     write(*,*) "written !"
!     read(*,*)
    
    call DMUMPS(imp%sys)
    imp_back%RHS(:) = imp%sys%RHS(:)
    
    do i = 1, imp%sys%NZ
        row = imp%sys%IRN(i)
        col = imp%sys%JCN(i)
        imp_back%ANZ(i) = -imp_back%RHS(row) * imp%RHS(col)
    end do
!     imp_back%RHS(:) = 0.0_rp
   
end subroutine solve_system_mumps_back
