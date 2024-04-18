!======================================================================================================================!
!
!                    DassFlow1D Version 2.1
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
!> @file preissmann_time_step.f90
!> @brief This file includes the subroutine for solving the Shallow-Water Equations using Preissmann scheme with LPI 
!> and double sweep methods.


!> \brief Solve the Shallow-Water Equations using Preissmann scheme with LPI and double sweep methods.
subroutine init_implicit_reduced_system(mdl, msh, imp)
    use m_common
    use m_mesh
    use m_sw_mono, only: ImplicitMatrix, Model, Unknowns
    implicit none
#ifndef CPP_ADJ
    include "mpif.h"
    include 'dmumps_struc.h'
#endif
    
    !=ARGUMENTS========================================================================================================!
    !> @param msh Model
    TYPE(Model), intent(in) :: mdl
    !> @param msh Mesh of the model
    TYPE(Mesh), intent(inout) :: msh
    !> @param[in] imp Implicit matrix
    TYPE(ImplicitMatrix), intent(inout) :: imp
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of column
    integer(ip) :: col
    ! Iterator
    integer(ip) :: i
    ! Index of boundary condition
    integer(ip) :: ibc
    ! Error code
    integer(ip) :: ierr
    ! Index of segment
    integer(ip) :: iseg
    ! Iterator on upstream segments
    integer(ip) :: ius_seg
    ! Size of the reduced system
    integer(ip) :: n
    ! Number of non-zero entries
    integer(ip) :: nnz
    ! Offset in non zero entries
    integer(ip) :: off
    ! Index of row
    integer(ip) :: row
    ! Index of upstream segment
    integer(ip) :: us_seg
    
    ! Initialise MUMPS
    call MPI_INIT(ierr)
    imp%sys%COMM = MPI_COMM_WORLD
    imp%sys%JOB = -1
    imp%sys%SYM = 0
    imp%sys%PAR = 1
    imp%sys%ICNTL(14) = 40
    call DMUMPS(imp%sys)
    imp%sys%ICNTL(4) = 1
    imp%sys%ICNTL(1) = 0
    imp%sys%ICNTL(2) = 0
    imp%sys%ICNTL(3) = 0
    
    ! Compute non-zero entries
    n = msh%nseg * 4
    nnz = 0
    do iseg = 1, msh%nseg
    
        ! Count non zero entries for upstream boundary condition
        if (msh%seg(iseg)%us_bc > 0) then
            ibc = msh%seg(iseg)%us_bc
            if (mdl%bc(ibc)%id == "discharge") then
                nnz = nnz + 1
            end if
        else
            nnz = nnz + size(msh%seg(iseg)%us_seg) + 1
        end if
    
        ! Count non zero entries for condensed transfers relations
        nnz = nnz + 8
        
        ! Count non zero entries for downstream boundary condition
        if (msh%seg(iseg)%ds_bc > 0) then
            ibc = msh%seg(iseg)%ds_bc
            if (mdl%bc(ibc)%id == "elevation") then
                nnz = nnz + 1
            else if (mdl%bc(ibc)%id == "normal_depth") then
                nnz = nnz + 2
            else if (mdl%bc(ibc)%id == "rating_curve") then
                nnz = nnz + 2
            end if
        else
            nnz = nnz + 2
        end if
        
    end do
    
    ! Allocate MUMPS arrays
    imp%sizes(1) = n
    imp%sizes(2) = nnz
    imp%sys%NZ = nnz
    imp%sys%N = n
    allocate(imp%ANZ(nnz))
    allocate(imp%sys%A(nnz))
    allocate(imp%sys%IRN(nnz))
    allocate(imp%sys%JCN(nnz))
    allocate(imp%RHS(n))
!     imp%sys%RHS => imp%RHS
    allocate(imp%sys%RHS(n))
    
    ! Compute rows and columns indices
    off = 1
    row = 1
    do iseg = 1, msh%nseg
    
        ! Upstream boundary condition
        if (msh%seg(iseg)%us_bc > 0) then
        
            ibc = msh%seg(iseg)%us_bc
            
            if (mdl%bc(ibc)%id == "discharge") then
            
                col = 4 * (iseg - 1) + 1
                imp%sys%IRN(off) = row
                imp%sys%JCN(off) = col
                off = off + 1
                
            end if
            
        else
        
            do ius_seg = 1, size(msh%seg(iseg)%us_seg)
                us_seg = msh%seg(iseg)%us_seg(ius_seg)
                col = 4 * (us_seg - 1) + 3
                imp%sys%IRN(off) = row
                imp%sys%JCN(off) = col
                off = off + 1
            end do
            col = 4 * (iseg - 1) + 1
            imp%sys%IRN(off) = row
            imp%sys%JCN(off) = col
            off = off + 1
            
        end if
        row = row + 1
        
        ! Condensend transfer relations
        col = 4 * (iseg - 1) + 1
        do i = 1, 4
            imp%sys%IRN(off) = row
            imp%sys%JCN(off) = col
            off = off + 1
            col = col + 1
        end do
        row = row + 1
        col = 4 * (iseg - 1) + 1
        do i = 1, 4
            imp%sys%IRN(off) = row
            imp%sys%JCN(off) = col
            off = off + 1
            col = col + 1
        end do
        row = row + 1
        
        ! Downstream boundary condition
        if (msh%seg(iseg)%ds_bc > 0) then
            ibc = msh%seg(iseg)%ds_bc
            if (mdl%bc(ibc)%id == "elevation") then
                col = 4 * (iseg - 1) + 4
                imp%sys%IRN(off) = row
                imp%sys%JCN(off) = col
                off = off + 1
            else if (mdl%bc(ibc)%id == "normal_depth") then
                col = 4 * (iseg - 1) + 3
                imp%sys%IRN(off) = row
                imp%sys%JCN(off) = col
                imp%sys%IRN(off+1) = row
                imp%sys%JCN(off+1) = col+1
                off = off + 2
            else if (mdl%bc(ibc)%id == "rating_curve") then
                col = 4 * (iseg - 1) + 3
                imp%sys%IRN(off) = row
                imp%sys%JCN(off) = col
                imp%sys%IRN(off+1) = row
                imp%sys%JCN(off+1) = col+1
                off = off + 2
            end if
        else
            col = 4 * (iseg - 1) + 4
            imp%sys%IRN(off) = row
            imp%sys%JCN(off) = col
            col = 4 * (msh%seg(iseg)%ds_seg - 1) + 2
            imp%sys%IRN(off+1) = row
            imp%sys%JCN(off+1) = col
            off = off + 2
        end if
        row = row + 1

    end do
    
!     imp%ANZ => imp%sys%A
!     imp%RHS => imp%sys%RHS
   
end subroutine init_implicit_reduced_system

