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
!> \file write_mesh.f90
!> \brief This file includes the routine for writing mesh
subroutine write_mesh(filename, msh)
    use m_mesh
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Mesh
    character(len=*), intent(in) :: filename
    !> Mesh
    type(Mesh), intent(in) :: msh
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of level
    integer(ip) :: ilevel
    ! Index of segment
    integer(ip) :: iseg
    ! Number of interior cross-sections
    integer(ip) :: ncs_int
    ! New index (for interior cross-sections)
    integer(ip) :: new_ics
    ! Number of levels
    integer(ip) :: nlevels
    ! Overbanks flag
    integer(ip) :: ob_flag
    ! Indices of interior cross-sections
    integer(ip), dimension(:), allocatable :: ics_int
    ! Error message
    character(len=lchar) :: err_msg
    
    ! Compute indices of interior cells
    allocate(ics_int(msh%ncs))
    if (msh%has_ghost_cells) then
        
        ! Set indices of all interior cells to 1 and ghost to 0
        do iseg = 1, msh%nseg
            ics_int(msh%seg(iseg)%first_cs-2:msh%seg(iseg)%first_cs-1) = 0
            ics_int(msh%seg(iseg)%first_cs:msh%seg(iseg)%last_cs) = 1
            ics_int(msh%seg(iseg)%last_cs+1:msh%seg(iseg)%last_cs+2) = 0
        end do
        
        ! Loop on all cross-sections to compute interior indices
        new_ics = 0
        do ics = 1, msh%ncs
            if (ics_int(ics) == 1) then
                new_ics = new_ics + 1
                ics_int(ics) = new_ics
            end if
        end do
        ncs_int = new_ics
        
    else
    
        do ics = 1, msh%ncs
            ics_int(ics) = ics
        end do
        ncs_int = msh%ncs
        
    end if
   
    ! Open file
    open(10, file=filename)
    
    ! Write header line
    write(10, '(A)', err=100) "DASSFLOW-1D MESH V1.3"
    write(10, '(A)', err=100) "# Generated using Python API"
    
    ! Write number of cross-sections and segments
    write(10, '(I4, I4, A)', err=100) ncs_int, msh%nseg, " linear none"
    
    ! Skip cross-section header line
    write(10, '(A)', err=100) "# CROSS-SECTIONS"

    do ics = 1, msh%ncs
    
        ! Skip ghost cells
        if (ics_int(ics) == 0) cycle
    
        ! Write global definition
        nlevels = size(msh%cs(ics)%level_heights)
        write(10, *, err=100) ics_int(ics), msh%cs(ics)%x, msh%cs(ics)%coord%x, msh%cs(ics)%coord%y, msh%cs(ics)%bathy, nlevels
        
        ! Write levels definition
        do ilevel = 1, nlevels
            ob_flag = 0
            if (msh%cs(ics)%ob_levels(2) <= ilevel .and. msh%cs(ics)%ob_levels(2) > 0) then
                ob_flag = ob_flag + 1
            end if
            if (msh%cs(ics)%ob_levels(1) <= ilevel .and. msh%cs(ics)%ob_levels(1) > 0) then
                ob_flag = ob_flag + 10
            end if
            write(10, *, err=100) msh%cs(ics)%level_heights(ilevel), msh%cs(ics)%level_widths(ilevel), 1.0, ob_flag
        end do
        
    end do
    
    ! Write segments
    write(10, '(A)', err=100) "# SEGMENTS"
    do iseg = 1, msh%nseg
        
        ! Write global definition
        write(10, *, err=100) iseg, ics_int(msh%seg(iseg)%first_cs), ics_int(msh%seg(iseg)%last_cs), &
                              size(msh%seg(iseg)%us_seg)
            
        ! Write indices of upstream and downstream segments
        write(10, *, err=100) msh%seg(iseg)%us_seg(:), msh%seg(iseg)%ds_seg
        
    end do
        
    close(10)
    
    return
  
100 write(err_msg, '(4A)') "error writing in file ", filename
    close(10)
    call f90wrap_abort(err_msg)
  
end subroutine
