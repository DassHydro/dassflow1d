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
!> @file curvilinear_abscissae.f90
!> @brief This file includes the routine that solves the Steady State Shallow Water Equation using standard step method

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  TODO
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> @brief Solve the Steady State Shallow Water Equation using standard step method
subroutine curvilinear_abscissae_from_coords(msh)
    use m_mesh, only: Mesh
    use m_mesh, only: ip, rp
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Mesh) :: msh
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of downstream segment
    integer(ip) :: ds_seg
    ! Index
    integer(ip) :: i
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of segment
    integer(ip) :: iseg
    ! Index of segment
    integer(ip) :: seg_stack_size
    ! Stack of segments
    integer(ip), dimension(msh%nseg) :: seg_stack
    ! Spacing in x direction
    real(rp) :: dx
    ! Spacing in y direction
    real(rp) :: dy

    ! Retrieve downstream segment
    ds_seg = 0
    do iseg = 1, msh%nseg
    
        if (msh%seg(iseg)%ds_seg <= 0) then
        
            ! Check that a previous downstream segment has not been found
            if (ds_seg /= 0) then
                call abort_solver("Multiple downstream segments in mesh")
            end if
            
            ! Set downstream segment index
            ds_seg = iseg
            
        end if
    
    end do
    
    if (ds_seg == 0) then
        call abort_solver("No downstream segment found in mesh")
    end if
    
    ! Initialise stack
    seg_stack(:) = 0
    seg_stack(1) = ds_seg
    seg_stack_size = 1
    
    ! Main loop
    do while(seg_stack_size > 0)
    
        ! Retrieve first segment in stack
        iseg = seg_stack(1)
        
        ! Retrieve downstream curvilinear abscissa
        ics = msh%seg(iseg)%last_cs
        if (msh%seg(iseg)%ds_seg <= 0) then
        
            msh%cs(ics)%x = 0.0
            
        else
        
            ds_seg = msh%seg(iseg)%ds_seg
            msh%cs(ics)%x = msh%cs(msh%seg(ds_seg)%first_cs)%x
            
        end if
        
        
        ! Compute curvilinear abscissae from downstream to upstream
        do ics = msh%seg(iseg)%last_cs - 1, msh%seg(iseg)%first_cs, -1
        
        
            dx = msh%cs(ics+1)%coord%x - msh%cs(ics)%coord%x
            dy = msh%cs(ics+1)%coord%y - msh%cs(ics)%coord%y
            msh%cs(ics)%x = msh%cs(ics+1)%x + sqrt(dx**2 + dy**2)
            msh%cs(ics+1)%deltademi = sqrt(dx**2 + dy**2)
        
        end do
        
        ! Compute curvilinear abscissa for ghost cross-sections
        msh%cs(msh%seg(iseg)%first_cs-1)%x =msh%cs(msh%seg(iseg)%first_cs)%x
        msh%cs(msh%seg(iseg)%first_cs-2)%x =msh%cs(msh%seg(iseg)%first_cs)%x
        msh%cs(msh%seg(iseg)%last_cs+1)%x =msh%cs(msh%seg(iseg)%last_cs)%x
        msh%cs(msh%seg(iseg)%last_cs+2)%x =msh%cs(msh%seg(iseg)%last_cs)%x
        
        ! Pop current segment from stack
        if (seg_stack_size > 1) then
          seg_stack(1:seg_stack_size-1) = seg_stack(2:seg_stack_size)
          seg_stack(seg_stack_size) = 0
        end if
        seg_stack_size = seg_stack_size - 1
        
        ! Append stack with upstream segments
        if (allocated(msh%seg(iseg)%us_seg)) then
          
            do i = 1, size(msh%seg(iseg)%us_seg)
                if (msh%seg(iseg)%us_seg(i) > 0) then
                    seg_stack_size = seg_stack_size + 1
                    seg_stack(seg_stack_size) = msh%seg(iseg)%us_seg(i)
                end if
            end do
            
        end if
        
    end do
    
end subroutine


subroutine finalise_curvilinear_abscissae(msh)
    use m_mesh, only: Mesh
    use m_mesh, only: ip, rp
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Mesh) :: msh
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of segment
    integer(ip) :: iseg
    
    do iseg = 1, msh%nseg
    
        ! Compute deltademi
        do ics = msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs - 1
            msh%cs(ics+1)%deltademi = abs(msh%cs(ics+1)%x - msh%cs(ics)%x)
        end do
        
        ! Compute curvilinear abscissa for ghost cross-sections
        msh%cs(msh%seg(iseg)%first_cs-1)%x =msh%cs(msh%seg(iseg)%first_cs)%x
        msh%cs(msh%seg(iseg)%first_cs-2)%x =msh%cs(msh%seg(iseg)%first_cs)%x
        msh%cs(msh%seg(iseg)%last_cs+1)%x =msh%cs(msh%seg(iseg)%last_cs)%x
        msh%cs(msh%seg(iseg)%last_cs+2)%x =msh%cs(msh%seg(iseg)%last_cs)%x
    
    end do
    
end subroutine
