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
!> \file read_mesh.f90
!> \brief This file includes the routine for reading mesh

! #define DEBUG_L1
! #define DEBUG_L2
!#define RECORD_INTERPOLATIONS

!> Resample a mesh
subroutine resample_mesh(msh, dxs)
    use m_common
!     use m_mesh, only: Mesh
!     use m_mesh, only: mesh_initialise, crosssection_copy
!     use m_mesh, only: CrossSection
    use m_mesh
    implicit none
    interface
        subroutine interpolate_between_crosssections(up_cs, dn_cs, n, interp_css)
            import CrossSection, ip
            type(CrossSection), intent(in) :: up_cs
            type(CrossSection), intent(in) :: dn_cs
            integer(ip), intent(in) :: n
            type(CrossSection), dimension(:), intent(inout) :: interp_css
        end subroutine
    end interface
    !=ARGUMENTS========================================================================================================!
    !> Mesh
    type(Mesh), intent(inout) :: msh
    !> Mesh
    real(rp) :: dxs
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of first cross-section
    integer(ip) :: first_cs
    ! Index
    integer(ip) :: i
    ! Index of base cross-section
    integer(ip) :: base_ics
    ! Index of cross-section
    integer(ip) :: ics
!     ! Index of level
!     integer(ip) :: ilevel
!     ! Index of line in mesh file
!     integer(ip) :: iline
    ! Index of segment
    integer(ip) :: iseg
!     ! Index of upstream segment
!     integer(ip) :: ius_seg
!     ! Index of last cross-section
!     integer(ip) :: last_cs
!     ! String for reading whole line in mesh file
!     character(len=256) :: line
!     ! Shape model
!     character(len=lchar) :: shape_model
!     ! Number of boundary conditions
!     integer(ip) :: nbc
    ! Number of cross-sections
    integer(ip) :: ncs
    ! Number of sections after interpolation between two cross-sections (included)
    integer(ip) :: ninterp
!     ! Number of levels
!     integer(ip) :: nlevels
!     ! Number of segments
!     integer(ip) :: nseg
!     ! Number of upstream segments indices
!     integer(ip) :: nus_seg
!     ! Overbanks flag
!     integer(ip) :: ob_flag
!     ! Projection
!     character(len=lchar) :: proj
!     ! Temporary value
!     real(rp) :: tmp
!     ! Interpolation coefficient
!     real(rp) :: alpha
    ! Recorded base cross-sections
    type(CrossSection), dimension(:), allocatable :: base_cs
!     
!     character(len=lchar) :: err_msg
    
    ! Add ghost_cells if necessary
    if (.not. msh%has_ghost_cells) call add_ghost_cells(msh)
    
    ! Store base cross-sections
    allocate(base_cs(msh%ncs)) 
    do ics = 1, msh%ncs
        call crosssection_copy(msh%cs(ics), base_cs(ics))
    end do
    
    ! Compute number of cross-sections after resampling
    ncs = 0
    do iseg = 1, msh%nseg
    
#ifdef DEBUG_L1
        print *, "iseg:", iseg, msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
#endif
    
        ! Add upstream ghost cross-sections
        ncs = ncs + 2
        
        ! Add interior cross-sections
        do ics = msh%seg(iseg)%first_cs + 1, msh%seg(iseg)%last_cs
#ifdef DEBUG_L1
            print *, ics, abs(base_cs(ics)%x - base_cs(ics-1)%x), &
                     max(1, nint(abs(base_cs(ics)%x - base_cs(ics-1)%x) / dxs) - 1), &
                     ncs, ncs + max(1, nint(abs(base_cs(ics)%x - base_cs(ics-1)%x) / dxs) - 1)
!                      
!                      nint(abs(base_cs(ics)%x - base_cs(ics-1)%x) / dxs)
#endif
            ncs = ncs + max(1, nint(abs(base_cs(ics)%x - base_cs(ics-1)%x) / dxs) - 1)
        end do
        ncs = ncs + 1
        
        ! Add downstream ghost cross-sections
        ncs = ncs + 2
    
    end do
    
#ifdef DEBUG_L1
    print *, "New number of cross-sections:", ncs
#endif
    
    ! Reallocate mesh with actual number of cross-sections (real + ghost)
    deallocate(msh%cs)
    allocate(msh%cs(ncs))
    msh%ncs = ncs
    
    ! Compute actual cross-sections after resampling
    ics = 1
    do iseg = 1, msh%nseg
    
        ! Add upstream ghost cross-sections
#ifdef DEBUG_L2
        print *, "=> append upstream ghost:", ics
#endif
        call crosssection_copy(base_cs(msh%seg(iseg)%first_cs-2), msh%cs(ics))
#ifdef DEBUG_L2
        print *, "=> append upstream ghost:", ics+1
#endif
        call crosssection_copy(base_cs(msh%seg(iseg)%first_cs-1), msh%cs(ics+1))
        ics = ics + 2
        first_cs = ics
        
        ! Add interior cross-sections
#ifdef DEBUG_L2
        print *, "=> append most upstream cross-section:", ics
#endif
        call crosssection_copy(base_cs(msh%seg(iseg)%first_cs), msh%cs(ics))
#ifdef DEBUG_L2
        print '(2I4, 2(1X,F9.3))', ics, msh%seg(iseg)%first_cs, msh%cs(ics)%x, msh%cs(ics)%bathy
#endif
        ics = ics + 1
        do base_ics = msh%seg(iseg)%first_cs + 1, msh%seg(iseg)%last_cs
            ninterp = nint(abs(base_cs(base_ics)%x - base_cs(base_ics-1)%x) / dxs)
            if (ninterp < 2) ninterp = 2
#ifdef DEBUG_L1
            print *, base_ics, ics, ninterp
#endif
#ifdef DEBUG_L2
            print *, "II", base_ics, abs(base_cs(base_ics)%x - base_cs(base_ics-1)%x), ninterp
            print *, "=> interpolated ics:", ics, ics+ninterp-3
#endif
            if (ninterp > 2) then
                call interpolate_between_crosssections(base_cs(base_ics-1), base_cs(base_ics), ninterp, &
                                                       msh%cs(ics:ics+ninterp-3))
            end if
#ifdef DEBUG_L2
            do i = 0, ninterp-3
                print '(I4,A,2(1X,F9.3))', ics+i, " ###", msh%cs(ics+i)%x, msh%cs(ics+i)%bathy
            end do
#endif
            ics = ics + ninterp - 2
#ifdef DEBUG_L2
            print *, "=> append existing downstream section:", ics, base_ics
#endif
            call crosssection_copy(base_cs(base_ics), msh%cs(ics))
#ifdef DEBUG_L2
            print '(2I4, 2(1X,F9.3))', ics, base_ics, msh%cs(ics)%x, msh%cs(ics)%bathy
#endif
            ics = ics + 1

        end do
        
        ! Add downstream ghost cross-sections
        call crosssection_copy(base_cs(msh%seg(iseg)%last_cs+1), msh%cs(ics))
        call crosssection_copy(base_cs(msh%seg(iseg)%last_cs+2), msh%cs(ics+1))
        ics = ics + 2
        
        msh%seg(iseg)%first_cs = first_cs
        msh%seg(iseg)%last_cs = ics-3

    
    end do
    
    call finalise_curvilinear_abscissae(msh)

    ! Compute slopes
    call bathy_slopes(msh)

    ! Check mesh
    call mesh_check(msh)
    
!     return
!   
! 100 write(buffer, '(I32)') iline
!     write(err_msg, '(4A)') "error in file ", filename, ", at line ", trim(adjustl(buffer))
!     close(10)
!     call f90wrap_abort(err_msg)
! 101 write(buffer, '(I32)') iline
!     write(err_msg, '(5A)') "error in file ", filename, ", at line ", trim(adjustl(buffer)), &
!                            " (indices of upstream segments must be lower than index of current segment)"
!     close(10)
!     call f90wrap_abort(err_msg)
! 102 write(buffer, '(I32)') iline
!     write(err_msg, '(5A)') "error in file ", filename, ", at line ", trim(adjustl(buffer)), &
!                            " (index of downstream segment must be greater than index of current segment)"
!     close(10)
!     call f90wrap_abort(err_msg)
  
end subroutine


! !> Resample a mesh
! subroutine interpolate_between_crosssections(up_cs, dn_cs, n, interp_css)
!     use m_mesh
!     use m_numeric
!     implicit none
!     !=ARGUMENTS========================================================================================================!
!     !> Upstream cross-section
!     type(CrossSection), intent(in) :: up_cs
!     !> Downstream cross-section
!     type(CrossSection), intent(in) :: dn_cs
!     !> Number of total cross-sections (existing + interpolated)
!     integer(ip), intent(in) :: n
!     !> Interpolated cross-sections
!     type(CrossSection), dimension(:), intent(inout) :: interp_css
!     !=LOCAL VARIABLES==================================================================================================!
!     ! Index
!     integer(ip) :: i
!     ! Index of cross-section
!     integer(ip) :: ics
!     ! Index of interpolated cross-section
!     integer(ip) :: iinterp
!     ! Index of segment
!     integer(ip) :: iseg
!     ! Indices of downstream master levels
!     integer(ip), dimension(3) :: dn_ilevel
!     ! Number of levels between downstream master levels
!     integer(ip), dimension(3) :: dn_nlevels
!     ! Indices of downstream overbank levels
!     integer(ip), dimension(2) :: dn_oblevels
!     ! Indices of master levels for interpolated cross-section
!     integer(ip), dimension(3) :: ilevel
!     ! Number of levels between master levels for interpolated cross-section
!     integer(ip), dimension(3) :: nlevels
!     ! Indices of overbank levels
!     integer(ip), dimension(2) :: ob_levels
!     ! Indices of upstream master levels
!     integer(ip), dimension(3) :: up_ilevel
!     ! Number of levels between upstream master levels
!     integer(ip), dimension(3) :: up_nlevels
!     ! Indices of upstream overbank levels
!     integer(ip), dimension(2) :: up_oblevels
! !     ! Index of last cross-section
! !     integer(ip) :: last_cs
! !     ! String for reading whole line in mesh file
! !     character(len=256) :: line
! !     ! Shape model
! !     character(len=lchar) :: shape_model
! !     ! Number of boundary conditions
! !     integer(ip) :: nbc
! !     ! Number of cross-sections
! !     integer(ip) :: ncs
! !     ! Number of levels
! !     integer(ip) :: nlevels
! !     ! Number of segments
! !     integer(ip) :: nseg
! !     ! Number of upstream segments indices
! !     integer(ip) :: nus_seg
! !     ! Overbanks flag
! !     integer(ip) :: ob_flag
! !     ! Projection
! !     character(len=lchar) :: proj
!     ! Interpolation coefficient
!     real(rp) :: alpha
!     ! Resampled downstream heights
!     real(rp), dimension(:), allocatable :: dn_heights
!     ! Resampled downstream widths
!     real(rp), dimension(:), allocatable :: dn_widths
!     ! Level heights for intepolated cross-section
!     real(rp), dimension(:), allocatable :: level_heights
!     ! Level widths for intepolated cross-section
!     real(rp), dimension(:), allocatable :: level_widths
!     ! Resampled upstream heights
!     real(rp), dimension(:), allocatable :: up_heights
!     ! Resampled upstream widths
!     real(rp), dimension(:), allocatable :: up_widths
! !     ! Mesh version
! !     real(rp) :: version
! !     ! Temporary cross-sections
! !     type(CrossSection), dimension(:), allocatable :: cs
! !     
! !     character(len=lchar) :: err_msg
!     
!     
!     ! Compute number and indices of master levels at upstream cross-section
!     if (up_cs%ob_levels(1) > 0 .and. up_cs%ob_levels(2) > 0) then
!         up_ilevel(1) = min(up_cs%ob_levels(1), up_cs%ob_levels(2))
!         up_ilevel(2) = max(up_cs%ob_levels(1), up_cs%ob_levels(2))
!         up_ilevel(3) = size(up_cs%level_heights)
!         if (up_cs%ob_levels(1) > up_cs%ob_levels(2)) then
!             up_oblevels(1) = up_ilevel(2)
!             up_oblevels(2) = up_ilevel(1)
!         else
!             up_oblevels(1) = up_ilevel(1)
!             up_oblevels(2) = up_ilevel(2)
!         end if
!     else if (up_cs%ob_levels(1) > 0) then
!         up_ilevel(1) = up_cs%ob_levels(1)
!         up_ilevel(2) = size(up_cs%level_heights)
!         up_ilevel(3) = up_ilevel(2)
!         up_oblevels(1) = up_ilevel(1)
!         up_oblevels(2) = up_ilevel(2)
!     else if (up_cs%ob_levels(2) > 0) then
!         up_ilevel(1) = up_cs%ob_levels(2)
!         up_ilevel(2) = size(up_cs%level_heights)
!         up_ilevel(3) = up_ilevel(2)
!         up_oblevels(1) = up_ilevel(2)
!         up_oblevels(2) = up_ilevel(1)
!     else
!         up_ilevel(1) = size(up_cs%level_heights)
!         up_ilevel(2) = up_ilevel(1)
!         up_ilevel(3) = up_ilevel(2)
!         up_oblevels(1) = up_ilevel(1)
!         up_oblevels(2) = up_ilevel(1)
!     end if
!     up_nlevels(1) = up_ilevel(1)
!     up_nlevels(2) = up_ilevel(2) - up_ilevel(1)
!     up_nlevels(3) = up_ilevel(3) - up_ilevel(2)
!     
!     ! Compute number and indices of master levels at downstream cross-section
!     if (dn_cs%ob_levels(1) > 0 .and. dn_cs%ob_levels(2) > 0) then
!         dn_ilevel(1) = min(dn_cs%ob_levels(1), dn_cs%ob_levels(2))
!         dn_ilevel(2) = max(dn_cs%ob_levels(1), dn_cs%ob_levels(2))
!         dn_ilevel(3) = size(dn_cs%level_heights)
!         if (dn_cs%ob_levels(1) > dn_cs%ob_levels(2)) then
!             dn_oblevels(1) = dn_ilevel(2)
!             dn_oblevels(2) = dn_ilevel(1)
!         else
!             dn_oblevels(1) = dn_ilevel(1)
!             dn_oblevels(2) = dn_ilevel(2)
!         end if
!     else if (dn_cs%ob_levels(1) > 0) then
!         dn_ilevel(1) = dn_cs%ob_levels(1)
!         dn_ilevel(2) = size(dn_cs%level_heights)
!         dn_ilevel(3) = dn_ilevel(2)
!         dn_oblevels(1) = dn_ilevel(1)
!         dn_oblevels(2) = dn_ilevel(2)
!     else if (dn_cs%ob_levels(2) > 0) then
!         dn_ilevel(1) = dn_cs%ob_levels(2)
!         dn_ilevel(2) = size(dn_cs%level_heights)
!         dn_ilevel(3) = dn_ilevel(2)
!         dn_oblevels(1) = dn_ilevel(2)
!         dn_oblevels(2) = dn_ilevel(1)
!     else
!         dn_ilevel(1) = size(dn_cs%level_heights)
!         dn_ilevel(2) = dn_ilevel(1)
!         dn_ilevel(3) = dn_ilevel(2)
!         dn_oblevels(1) = dn_ilevel(1)
!         dn_oblevels(2) = dn_ilevel(1)
!     end if
!     dn_nlevels(1) = dn_ilevel(1)
!     dn_nlevels(2) = dn_ilevel(2) - dn_ilevel(1)
!     dn_nlevels(3) = dn_ilevel(3) - dn_ilevel(2)
!     
! #ifdef DEBUG_L1
!     print *, "LEVELS:", up_ilevel(1), dn_ilevel(1)
!     print *, "       ", up_ilevel(2), dn_ilevel(2)
!     print *, "       ", up_ilevel(3), dn_ilevel(3)
! #endif
!             
!     ! Interpolate between cross-sections
!     do iinterp = 1, n-2
!         
!         alpha = real(iinterp) / real(n-1)
! #ifdef DEBUG_L1
!         print *, "ALPHA:", iinterp, alpha
! #endif
!         
!         ! Compute number of levels
!         nlevels = nint(up_cs%levels + alpha * (dn_cs%levels - up_cs%levels))
!         
!         !Resample upstream levels and heights
!         allocate(up_heights(nlevels))
!         allocate(up_widths(nlevels(1) + nlevels(2) + nlevels(3)))
!         
!         
! #ifdef DEBUG_L1
!         print *, "=>LEVELS:", nlevels(1)
!         print *, "         ", up_ilevel(2), dn_ilevel(2)
!         print *, "         ", up_ilevel(3), dn_ilevel(3)
! #endif
!             
!         ! Interpolate levels at upstream
!         allocate(up_heights(nlevels(1) + nlevels(2) + nlevels(3)))
!         allocate(up_widths(nlevels(1) + nlevels(2) + nlevels(3)))
!         up_heights(1) = up_cs%level_heights(1)
!         up_heights(nlevels(1)) = up_cs%level_heights(up_ilevel(1))
!         do i = 2, nlevels(1)-1
!             up_heights(i) = up_heights(1) + real(i-1)/real(nlevels(1)-1) * (up_heights(ilevel(1)) - up_heights(1))
!             up_widths(i) = linear_interp(up_cs%level_heights, up_cs%level_widths, up_heights(i))
!         end do
!         
!         up_heights(ilevel(2)) = up_cs%level_heights(up_ilevel(2))
!         if (nlevels(2) > 2) then
!             do i = 2, nlevels(2)-1
!                 up_heights(ilevel(1)+i-1) = up_heights(ilevel(1)) + real(i-1)/real(nlevels(2)-1) * (up_heights(ilevel(2)) - up_heights(ilevel(1)))
!                 up_widths(ilevel(1)+i-1) = linear_interp(up_cs%level_heights, up_cs%level_widths, up_heights(ilevel(1)+i-1))
!             end do
!         end if
!         
!         up_heights(ilevel(3)) = up_cs%level_heights(up_ilevel(3))
!         if (nlevels(3) > 2) then
!             do i = 2, nlevels(3)-1
!                 up_heights(ilevel(2)+i-1) = up_heights(ilevel(2)) + real(i-1)/real(nlevels(3)-1) * (up_heights(ilevel(3)) - up_heights(ilevel(2)))
!                 up_widths(ilevel(2)+i-1) = linear_interp(up_cs%level_heights, up_cs%level_widths, up_heights(ilevel(2)+i-1))
!             end do
!         end if
!             
!         ! Interpolate levels at downstream
!         allocate(dn_heights(nlevels(1) + nlevels(2) + nlevels(3)))
!         allocate(dn_widths(nlevels(1) + nlevels(2) + nlevels(3)))
!         dn_heights(1) = dn_cs%level_heights(1)
!         dn_heights(nlevels(1)) = dn_cs%level_heights(dn_ilevel(1))
!         do i = 2, nlevels(1)-1
!             dn_heights(i) = dn_heights(1) + real(i-1)/real(nlevels(1)-1) * (dn_heights(ilevel(1)) - dn_heights(1))
!             dn_widths(i) = linear_interp(dn_cs%level_heights, dn_cs%level_widths, dn_heights(i))
!         end do
! #ifdef DEBUG_L2
!         print *, "dn_heights(1:n1):", dn_heights(1:nlevels(1))
!         read(*,*)
! #endif
! 
!         dn_heights(ilevel(2)) = dn_cs%level_heights(dn_ilevel(2))
!         if (nlevels(2) > 2) then
!             do i = 2, nlevels(2)-1
!                 dn_heights(ilevel(1)+i-1) = dn_heights(ilevel(1)) + real(i-1)/real(nlevels(2)-1) * (dn_heights(ilevel(2)) - dn_heights(ilevel(1)))
!                 dn_widths(ilevel(1)+i-1) = linear_interp(dn_cs%level_heights, dn_cs%level_widths, dn_heights(ilevel(1)+i-1))
!             end do
!         end if
! #ifdef DEBUG_L2
!         print *, "dn_heights(n1:n2):", dn_heights(nlevels(1):nlevels(1)+nlevels(2))
!         read(*,*)
! #endif
!         
!         dn_heights(ilevel(3)) = dn_cs%level_heights(dn_ilevel(3))
!         if (nlevels(3) > 2) then
!             do i = 2, nlevels(3)-1
!                 dn_heights(ilevel(2)+i-1) = dn_heights(ilevel(2)) + real(i-1)/real(nlevels(3)-1) * (dn_heights(ilevel(3)) - dn_heights(ilevel(2)))
!                 dn_widths(ilevel(2)+i-1) = linear_interp(dn_cs%level_heights, dn_cs%level_widths, dn_heights(ilevel(2)+i-1))
!             end do
!         end if
! #ifdef DEBUG_L2
!         print *, "dn_heights(n2:n3):", dn_heights(nlevels(1)+nlevels(2):nlevels(1)+nlevels(2)+nlevels(3))
!         read(*,*)
!         print *, "dn_heights:", dn_heights
!         read(*,*)
! #endif
!         
!         ! Allocate levels at interpolated cs
! !         interp_css(iinterp)%nlevels = nlevels(1) + nlevels(2) + nlevels(3)
!         allocate(level_heights(nlevels(1) + nlevels(2) + nlevels(3)))
!         allocate(level_widths(nlevels(1) + nlevels(2) + nlevels(3)))
! !         allocate(interp_css(iinterp)%y(nlevels(1) + nlevels(2) + nlevels(3)))
!         do i = 1, nlevels(1) + nlevels(2) + nlevels(3)
!             level_heights(i) = (1.0_rp - alpha) * up_heights(i) + alpha * dn_heights(i)
!             level_widths(i) = (1.0_rp - alpha) * up_widths(i) + alpha * dn_widths(i)
!         end do
!         
!         call set_levels(interp_css(iinterp), level_heights, level_widths)
!         
!         interp_css(iinterp)%coord%x = (1.0_rp - alpha) * up_cs%coord%x + alpha * dn_cs%coord%x
!         interp_css(iinterp)%coord%y = (1.0_rp - alpha) * up_cs%coord%y + alpha * dn_cs%coord%y
!         interp_css(iinterp)%x = (1.0_rp - alpha) * up_cs%x + alpha * dn_cs%x
!         interp_css(iinterp)%ibase = 0
!         interp_css(iinterp)%level = nint((1.0_rp - alpha) * up_cs%level + alpha * dn_cs%level)
!         interp_css(iinterp)%bathy = (1.0_rp - alpha) * up_cs%bathy + alpha * dn_cs%bathy
!         interp_css(iinterp)%ob_levels(:) = ob_levels(:)
! #ifdef FLOODPLAIN_MODEL
!         interp_css(iinterp)%alphaFP = (1.0_rp - alpha) * up_cs%alphaFP + alpha * dn_cs%alphaFP
! #endif
!         interp_css(iinterp)%slope = up_cs%slope
!         interp_css(iinterp)%strickler_type_code = up_cs%strickler_type_code
!         interp_css(iinterp)%strickler_params(:) = (1.0_rp - alpha) * up_cs%strickler_params(:) + &
!                                                   alpha * dn_cs%strickler_params(:)
!     
!         deallocate(up_heights)
!         deallocate(up_widths)
!         deallocate(level_heights)
!         deallocate(level_widths)
!         deallocate(dn_heights)
!         deallocate(dn_widths)
!         
!         call update_geometry(interp_css(iinterp))
!         
!     end do
! 
! end subroutine
! 
! 
! 








!> Resample a mesh
subroutine interpolate_between_crosssections(up_cs, dn_cs, n, interp_css)
    use m_mesh
    use m_numeric
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Upstream cross-section
    type(CrossSection), intent(in) :: up_cs
    !> Downstream cross-section
    type(CrossSection), intent(in) :: dn_cs
    !> Number of total cross-sections (existing + interpolated)
    integer(ip), intent(in) :: n
    !> Interpolated cross-sections
    type(CrossSection), dimension(:), intent(inout) :: interp_css
    !=LOCAL VARIABLES==================================================================================================!
    ! Index
    integer(ip) :: i
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of interpolated cross-section
    integer(ip) :: iinterp
    ! Index of segment
    integer(ip) :: iseg
    ! Indices of downstream master levels
    integer(ip), dimension(3) :: dn_ilevel
    ! Number of levels between downstream master levels
    integer(ip), dimension(3) :: dn_nlevels
    ! Indices of downstream overbank levels
    integer(ip), dimension(2) :: dn_oblevels
    ! Indices of master levels for interpolated cross-section
    integer(ip), dimension(3) :: ilevel
    ! Number of levels between master levels for interpolated cross-section
    integer(ip), dimension(3) :: nlevels
    ! Indices of overbank levels
    integer(ip), dimension(2) :: ob_levels
    ! Indices of upstream master levels
    integer(ip), dimension(3) :: up_ilevel
    ! Number of levels between upstream master levels
    integer(ip), dimension(3) :: up_nlevels
    ! Indices of upstream overbank levels
    integer(ip), dimension(2) :: up_oblevels
!     ! Index of last cross-section
!     integer(ip) :: last_cs
!     ! String for reading whole line in mesh file
!     character(len=256) :: line
!     ! Shape model
!     character(len=lchar) :: shape_model
!     ! Number of boundary conditions
!     integer(ip) :: nbc
!     ! Number of cross-sections
!     integer(ip) :: ncs
!     ! Number of levels
!     integer(ip) :: nlevels
!     ! Number of segments
!     integer(ip) :: nseg
!     ! Number of upstream segments indices
!     integer(ip) :: nus_seg
    ! Overbanks flag
    integer(ip) :: ob_flag
!     ! Projection
!     character(len=lchar) :: proj
    ! Interpolation coefficient
    real(rp) :: alpha
    ! Resampled downstream heights
    real(rp), dimension(:), allocatable :: dn_heights
    ! Resampled downstream widths
    real(rp), dimension(:), allocatable :: dn_widths
    ! Level heights for intepolated cross-section
    real(rp), dimension(:), allocatable :: level_heights
    ! Level widths for intepolated cross-section
    real(rp), dimension(:), allocatable :: level_widths
    ! Resampled upstream heights
    real(rp), dimension(:), allocatable :: up_heights
    ! Resampled upstream widths
    real(rp), dimension(:), allocatable :: up_widths
!     ! Mesh version
!     real(rp) :: version
!     ! Temporary cross-sections
!     type(CrossSection), dimension(:), allocatable :: cs
!     
!     character(len=lchar) :: err_msg

#ifdef RECORD_INTERPOLATIONS
    character(len=64) :: fname
    write(fname, '(2(A,I4.4),A)') "interp/cs", up_cs%ibase, "x", dn_cs%ibase, ".dat"
    open(24, file=fname)
    write(24, '(A)') "INDEX HEIGHT WIDTH OB"
    do i = 1, size(up_cs%level_heights)
        ob_flag = 0
        if (i >= up_cs%ob_levels(1)) then
            ob_flag = ob_flag + 10
        end if
        if (i >= up_cs%ob_levels(2)) then
            ob_flag = ob_flag + 1
        end if
        write(24, '(I4,2(1X,E15.6),1X,I2.2)') 0, up_cs%level_heights(i), up_cs%level_widths(i), ob_flag
    end do
#endif
    
    ! Compute number and indices of master levels at upstream cross-section
    if (up_cs%ob_levels(1) > 0 .and. up_cs%ob_levels(2) > 0) then
        up_ilevel(1) = min(up_cs%ob_levels(1), up_cs%ob_levels(2))
        up_ilevel(2) = max(up_cs%ob_levels(1), up_cs%ob_levels(2))
        up_ilevel(3) = size(up_cs%level_heights)
        if (up_cs%ob_levels(1) > up_cs%ob_levels(2)) then
            up_oblevels(1) = up_ilevel(2)
            up_oblevels(2) = up_ilevel(1)
        else
            up_oblevels(1) = up_ilevel(1)
            up_oblevels(2) = up_ilevel(2)
        end if
    else if (up_cs%ob_levels(1) > 0) then
        up_ilevel(1) = up_cs%ob_levels(1)
        up_ilevel(2) = size(up_cs%level_heights)
        up_ilevel(3) = up_ilevel(2)
        up_oblevels(1) = up_ilevel(1)
        up_oblevels(2) = up_ilevel(2)
    else if (up_cs%ob_levels(2) > 0) then
        up_ilevel(1) = up_cs%ob_levels(2)
        up_ilevel(2) = size(up_cs%level_heights)
        up_ilevel(3) = up_ilevel(2)
        up_oblevels(1) = up_ilevel(2)
        up_oblevels(2) = up_ilevel(1)
    else
        up_ilevel(1) = size(up_cs%level_heights)
        up_ilevel(2) = up_ilevel(1)
        up_ilevel(3) = up_ilevel(2)
        up_oblevels(1) = up_ilevel(1)
        up_oblevels(2) = up_ilevel(1)
    end if
    up_nlevels(1) = up_ilevel(1)
    up_nlevels(2) = up_ilevel(2) - up_ilevel(1)
    up_nlevels(3) = up_ilevel(3) - up_ilevel(2)
    
    ! Compute number and indices of master levels at downstream cross-section
    if (dn_cs%ob_levels(1) > 0 .and. dn_cs%ob_levels(2) > 0) then
        dn_ilevel(1) = min(dn_cs%ob_levels(1), dn_cs%ob_levels(2))
        dn_ilevel(2) = max(dn_cs%ob_levels(1), dn_cs%ob_levels(2))
        dn_ilevel(3) = size(dn_cs%level_heights)
        if (dn_cs%ob_levels(1) > dn_cs%ob_levels(2)) then
            dn_oblevels(1) = dn_ilevel(2)
            dn_oblevels(2) = dn_ilevel(1)
        else
            dn_oblevels(1) = dn_ilevel(1)
            dn_oblevels(2) = dn_ilevel(2)
        end if
    else if (dn_cs%ob_levels(1) > 0) then
        dn_ilevel(1) = dn_cs%ob_levels(1)
        dn_ilevel(2) = size(dn_cs%level_heights)
        dn_ilevel(3) = dn_ilevel(2)
        dn_oblevels(1) = dn_ilevel(1)
        dn_oblevels(2) = dn_ilevel(2)
    else if (dn_cs%ob_levels(2) > 0) then
        dn_ilevel(1) = dn_cs%ob_levels(2)
        dn_ilevel(2) = size(dn_cs%level_heights)
        dn_ilevel(3) = dn_ilevel(2)
        dn_oblevels(1) = dn_ilevel(2)
        dn_oblevels(2) = dn_ilevel(1)
    else
        dn_ilevel(1) = size(dn_cs%level_heights)
        dn_ilevel(2) = dn_ilevel(1)
        dn_ilevel(3) = dn_ilevel(2)
        dn_oblevels(1) = dn_ilevel(1)
        dn_oblevels(2) = dn_ilevel(1)
    end if
    dn_nlevels(1) = dn_ilevel(1)
    dn_nlevels(2) = dn_ilevel(2) - dn_ilevel(1)
    dn_nlevels(3) = dn_ilevel(3) - dn_ilevel(2)
    
#ifdef DEBUG_L2
    print *, "LEVELS:", up_ilevel(1), dn_ilevel(1)
    print *, "       ", up_ilevel(2), dn_ilevel(2)
    print *, "       ", up_ilevel(3), dn_ilevel(3)
#endif
            
    ! Interpolate between cross-sections
    do iinterp = 1, n-2
        
        alpha = real(iinterp) / real(n-1)
#ifdef DEBUG_L2
        print *, "ALPHA:", iinterp, alpha
#endif
        
        ! Compute number of levels
        nlevels(1) = nint(up_nlevels(1) + alpha * (dn_nlevels(1) - up_nlevels(1)))
        nlevels(2) = nint(up_nlevels(2) + alpha * (dn_nlevels(2) - up_nlevels(2)))
        nlevels(3) = nint(up_nlevels(3) + alpha * (dn_nlevels(3) - up_nlevels(3)))
        if (up_nlevels(1) == 2 .or. dn_nlevels(1) == 2) then
            nlevels(1) = 2
        end if
        ilevel(1) = nlevels(1)
        ilevel(2) = nlevels(1)+nlevels(2)
        ilevel(3) = nlevels(1)+nlevels(2)+nlevels(3)
        ob_levels(1) = nint(up_oblevels(1) + alpha * (dn_oblevels(1) - up_oblevels(1)))
        ob_levels(2) = nint(up_oblevels(2) + alpha * (dn_oblevels(2) - up_oblevels(2)))
#ifdef DEBUG_L2
        print *, "=>LEVELS:", nlevels(:)
        print *, "         ", up_ilevel(1), dn_ilevel(1)
        print *, "         ", up_ilevel(2), dn_ilevel(2)
        print *, "         ", up_ilevel(3), dn_ilevel(3)
#endif
            
        ! Interpolate levels at upstream
        allocate(up_heights(nlevels(1) + nlevels(2) + nlevels(3)))
        allocate(up_widths(nlevels(1) + nlevels(2) + nlevels(3)))
        up_heights(:) = 0.0_rp
        up_heights(1) = up_cs%level_heights(1)
        up_widths(1) = up_cs%level_widths(1)
        up_heights(nlevels(1)) = up_cs%level_heights(up_ilevel(1))
        up_widths(nlevels(1)) = up_cs%level_widths(up_ilevel(1))
        do i = 2, nlevels(1)-1
            up_heights(i) = up_heights(1) + real(i-1)/real(nlevels(1)-1) * (up_heights(ilevel(1)) - up_heights(1))
            up_widths(i) = linear_interp(up_cs%level_heights, up_cs%level_widths, up_heights(i))
        end do
        
        up_heights(ilevel(2)) = up_cs%level_heights(up_ilevel(2))
        up_widths(ilevel(2)) = up_cs%level_widths(up_ilevel(2))
        if (nlevels(2) > 1) then
            do i = 2, nlevels(2)
                up_heights(ilevel(1)+i-1) = up_heights(ilevel(1)) + real(i-1)/real(nlevels(2)) * (up_heights(ilevel(2)) - up_heights(ilevel(1)))
                up_widths(ilevel(1)+i-1) = linear_interp(up_cs%level_heights, up_cs%level_widths, up_heights(ilevel(1)+i-1))
            end do
        end if
        
        up_heights(ilevel(3)) = up_cs%level_heights(up_ilevel(3))
        up_widths(ilevel(3)) = up_cs%level_widths(up_ilevel(3))
        if (nlevels(3) > 1) then
            do i = 2, nlevels(3)
                up_heights(ilevel(2)+i-1) = up_heights(ilevel(2)) + real(i-1)/real(nlevels(3)) * (up_heights(ilevel(3)) - up_heights(ilevel(2)))
                up_widths(ilevel(2)+i-1) = linear_interp(up_cs%level_heights, up_cs%level_widths, up_heights(ilevel(2)+i-1))
            end do
        end if
            
        ! Interpolate levels at downstream
        allocate(dn_heights(nlevels(1) + nlevels(2) + nlevels(3)))
        allocate(dn_widths(nlevels(1) + nlevels(2) + nlevels(3)))
        dn_heights(:) = 0.0_rp
        dn_heights(1) = dn_cs%level_heights(1)
        dn_widths(1) = dn_cs%level_widths(1)
        dn_heights(nlevels(1)) = dn_cs%level_heights(dn_ilevel(1))
        dn_widths(nlevels(1)) = dn_cs%level_widths(dn_ilevel(1))
        do i = 2, nlevels(1)-1
            dn_heights(i) = dn_heights(1) + real(i-1)/real(nlevels(1)-1) * (dn_heights(ilevel(1)) - dn_heights(1))
            dn_widths(i) = linear_interp(dn_cs%level_heights, dn_cs%level_widths, dn_heights(i))
        end do
#ifdef DEBUG_L2
        print *, "dn_heights(1:n1):", dn_heights(1:nlevels(1))
        read(*,*)
#endif

        dn_heights(ilevel(2)) = dn_cs%level_heights(dn_ilevel(2))
        dn_widths(ilevel(2)) = dn_cs%level_widths(dn_ilevel(2))
        if (nlevels(2) > 1) then
            do i = 2, nlevels(2)
                dn_heights(ilevel(1)+i-1) = dn_heights(ilevel(1)) + real(i-1)/real(nlevels(2)) * (dn_heights(ilevel(2)) - dn_heights(ilevel(1)))
                dn_widths(ilevel(1)+i-1) = linear_interp(dn_cs%level_heights, dn_cs%level_widths, dn_heights(ilevel(1)+i-1))
            end do
        end if
#ifdef DEBUG_L2
        print *, "dn_heights(n1:n2):", dn_heights(nlevels(1):nlevels(1)+nlevels(2))
        read(*,*)
#endif
        
        dn_heights(ilevel(3)) = dn_cs%level_heights(dn_ilevel(3))
        dn_widths(ilevel(3)) = dn_cs%level_widths(dn_ilevel(3))
        if (nlevels(3) > 1) then
            do i = 2, nlevels(3)
                dn_heights(ilevel(2)+i-1) = dn_heights(ilevel(2)) + real(i-1)/real(nlevels(3)) * (dn_heights(ilevel(3)) - dn_heights(ilevel(2)))
                dn_widths(ilevel(2)+i-1) = linear_interp(dn_cs%level_heights, dn_cs%level_widths, dn_heights(ilevel(2)+i-1))
            end do
        end if
#ifdef DEBUG_L2
        print *, "dn_heights(n2:n3):", dn_heights(nlevels(1)+nlevels(2):nlevels(1)+nlevels(2)+nlevels(3))
        read(*,*)
        print *, "dn_heights:", dn_heights
        read(*,*)
#endif
        
        ! Allocate levels at interpolated cs
!         interp_css(iinterp)%nlevels = nlevels(1) + nlevels(2) + nlevels(3)
        allocate(level_heights(nlevels(1) + nlevels(2) + nlevels(3)))
        allocate(level_widths(nlevels(1) + nlevels(2) + nlevels(3)))
!         allocate(interp_css(iinterp)%y(nlevels(1) + nlevels(2) + nlevels(3)))
        do i = 1, nlevels(1) + nlevels(2) + nlevels(3)
            level_heights(i) = (1.0_rp - alpha) * up_heights(i) + alpha * dn_heights(i)
            level_widths(i) = (1.0_rp - alpha) * up_widths(i) + alpha * dn_widths(i)
#ifdef RECORD_INTERPOLATIONS
            ob_flag = 0
            if (i >= ob_levels(1)) then
                ob_flag = ob_flag + 10
            end if
            if (i >= ob_levels(2)) then
                ob_flag = ob_flag + 1
            end if
            write(24, '(I4,2(1X,E15.6),1X,I2.2)') iinterp, level_heights(i), level_widths(i), ob_flag
#endif
        end do
        
        call set_levels(interp_css(iinterp), level_heights, level_widths)
        
        interp_css(iinterp)%coord%x = (1.0_rp - alpha) * up_cs%coord%x + alpha * dn_cs%coord%x
        interp_css(iinterp)%coord%y = (1.0_rp - alpha) * up_cs%coord%y + alpha * dn_cs%coord%y
        interp_css(iinterp)%x = (1.0_rp - alpha) * up_cs%x + alpha * dn_cs%x
        interp_css(iinterp)%ibase = 0
        interp_css(iinterp)%level = nint((1.0_rp - alpha) * up_cs%level + alpha * dn_cs%level)
        interp_css(iinterp)%bathy = (1.0_rp - alpha) * up_cs%bathy + alpha * dn_cs%bathy
        interp_css(iinterp)%ob_levels(:) = ob_levels(:)
#ifdef FLOODPLAIN_MODEL
        interp_css(iinterp)%alphaFP = (1.0_rp - alpha) * up_cs%alphaFP + alpha * dn_cs%alphaFP
#endif
        interp_css(iinterp)%slope = up_cs%slope
        interp_css(iinterp)%strickler_type_code = up_cs%strickler_type_code
        interp_css(iinterp)%strickler_params(:) = (1.0_rp - alpha) * up_cs%strickler_params(:) + &
                                                  alpha * dn_cs%strickler_params(:)
    
        deallocate(up_heights)
        deallocate(up_widths)
        deallocate(level_heights)
        deallocate(level_widths)
        deallocate(dn_heights)
        deallocate(dn_widths)
        
        call update_geometry(interp_css(iinterp))
        
    end do
    
#ifdef RECORD_INTERPOLATIONS
    do i = 1, size(dn_cs%level_heights)
        ob_flag = 0
        if (i >= dn_cs%ob_levels(1)) then
            ob_flag = ob_flag + 10
        end if
        if (i >= dn_cs%ob_levels(2)) then
            ob_flag = ob_flag + 1
        end if
        write(24, '(I4,2(1X,E15.6),1X,I2.2)') -1, dn_cs%level_heights(i), dn_cs%level_widths(i), ob_flag
    end do
    close(24)
#endif

end subroutine
