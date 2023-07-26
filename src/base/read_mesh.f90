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

!> Compute the cumulated variable for all cross-sections of the msh
subroutine read_mesh(filename, msh)
    use m_mesh
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Mesh
    character(len=*), intent(in) :: filename
    !> Mesh
    type(Mesh), intent(out) :: msh
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of first cross-section
    integer(ip) :: first_cs
    ! Index
    integer(ip) :: i
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of level
    integer(ip) :: ilevel
    ! Index of line in mesh file
    integer(ip) :: iline
    ! Index of segment
    integer(ip) :: iseg
    ! Index of upstream segment
    integer(ip) :: ius_seg
    ! Index of last cross-section
    integer(ip) :: last_cs
    ! String for reading whole line in mesh file
    character(len=256) :: line
    ! Shape model
    character(len=lchar) :: shape_model
    ! Number of boundary conditions
    integer(ip) :: nbc
    ! Number of cross-sections
    integer(ip) :: ncs
    ! Number of levels
    integer(ip) :: nlevels
    ! Number of segments
    integer(ip) :: nseg
    ! Number of upstream segments indices
    integer(ip) :: nus_seg
    ! Overbanks flag
    integer(ip) :: ob_flag
    ! Projection
    character(len=lchar) :: proj
    ! Temporary value
    real(rp) :: tmp
    ! Mesh version
    real(rp) :: version
    ! Temporary cross-sections
    type(CrossSection), dimension(:), allocatable :: cs
    
    character(len=512) :: err_msg
   
    ! Open file
    open(10, file=filename)
    
    ! Read header line
    iline = 1
    read(10, '(A)', err=100, end=100) line
    if (line(1:18) == "DASSFLOW-1D MESH V") then
        ! Read mesh version
        read(line(19:len_trim(line)), *, err=100, end=100) version
!         print *, version
        
        ! Read header (comment) line
        iline = iline + 1
        read(10, '(A)', err=100, end=100) line
    else
        version = 0.0
    end if
    
    ! Read number of cross-sections and segments
    iline = iline + 1
    if (version < 1.0) then
        read(10, *, err=100, end=100) ncs, tmp
        nseg = 1
        shape_model = 'linear'
    else if (version < 1.1) then
        read(10, *, err=100, end=100) ncs, nseg
        shape_model = 'linear'
    else if (version < 1.3) then
        read(10, *, err=100, end=100) ncs, nseg, shape_model
    else
        read(10, *, err=100, end=100) ncs, nseg, shape_model, proj
    end if
    
    ! Skip cross-section header line
    iline = iline + 1
    read(10, '(A)', err=100, end=100) line
    
    ! Allocate and read temporary cross-sections
    allocate(cs(ncs))
    do i = 1, ncs
    
        ! Read global definition
        iline = iline + 1
        if (version < 1.3) then
            read(10, *, err=100, end=100) ics, cs(ics)%coord%x, cs(ics)%coord%y, cs(ics)%bathy, nlevels
            cs(ics)%ibase = ics
        else
            read(10, *, err=100, end=100) ics, cs(ics)%x, cs(ics)%coord%x, cs(ics)%coord%y, cs(ics)%bathy, nlevels
            cs(ics)%ibase = ics
        end if
        
        ! Read levels definition
        call crosssection_initialise(cs(ics), nlevels, shape_model)
        do ilevel = 1, nlevels
            iline = iline + 1
            if (version < 1.1) then
                read(10, *, err=100, end=100) cs(ics)%level_heights(ilevel), cs(ics)%level_widths(ilevel)
                ob_flag = 0
            else
                ! TODO: read center displacement (in tmp for now)
                read(10, *, err=100, end=100) cs(ics)%level_heights(ilevel), cs(ics)%level_widths(ilevel), tmp, ob_flag
            end if
            if (int(ob_flag/10) == 1 .and. cs(ics)%ob_levels(1) == 0) then
                cs(ics)%ob_levels(1) = ilevel
            end if
            if (mod(ob_flag,10) == 1 .and. cs(ics)%ob_levels(2) == 0) then
                cs(ics)%ob_levels(2) = ilevel
            end if
        end do
        
    end do
    
    ! Allocate segments
    allocate(msh%seg(nseg))
    msh%nseg = nseg
    
    ! Read or create segments
    if (version == 0.0) then
    
        ! Set the single segment for a DassFlow-1D v1.0 mesh file (segments not specified in the file)
        msh%seg(1)%first_cs = 1
        msh%seg(1)%last_cs = ncs
        allocate(msh%seg(1)%us_seg(1))
        nbc = 2
        msh%seg(1)%us_seg(1) = -1
        msh%seg(1)%us_bc = 1
        msh%seg(1)%ds_seg = -1
        msh%seg(1)%ds_bc = 2
        
    else
    
        ! Skip segments header line
        iline = iline + 1
        read(10, '(A)', err=100, end=100) line
        nbc = 0
        do i = 1, nseg
        
            ! Read global definition
            iline = iline + 1
            read(10, *, err=100, end=100) iseg, msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs, nus_seg
            
            ! Read indices of upstream and downstream segments
            allocate(msh%seg(iseg)%us_seg(nus_seg))
            iline = iline + 1
            read(10, *, err=100, end=100) msh%seg(iseg)%us_seg(:), msh%seg(iseg)%ds_seg
            
            ! Check that indices of upstream segments are lower than index of current segment
            do ius_seg = 1, nus_seg
                if (msh%seg(iseg)%us_seg(ius_seg) >= iseg) goto 101
            end do
            
            ! Check that index of downstream segment is greater than index of current segment
            if (msh%seg(iseg)%ds_seg > 0) then
                if (msh%seg(iseg)%ds_seg <= iseg) goto 102
            end if
            
            if (msh%seg(iseg)%us_seg(1) < 1) then
                nbc = nbc + 1
                msh%seg(iseg)%us_bc = nbc
            else
                msh%seg(iseg)%us_bc = 0
            end if
            if (msh%seg(iseg)%ds_seg < 1) then
                nbc = nbc + 1
                msh%seg(iseg)%ds_bc = nbc
            else
                msh%seg(iseg)%ds_bc = 0
            end if
            
        end do
        
        
    end if
    
    ! Close file
    close(10)
    
    ! Compute actual number of cross-sections (real + ghost)
    ncs = 0
    do iseg = 1, nseg
    
      ! Add number of upstream ghost cross-sections
      ncs = ncs + 2
      ! Add number of interior cross-sections
      ncs = ncs + msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs + 1
      ! Add number of downstream ghost cross-sections
      ncs = ncs + 2
      
    end do
    
    ! Allocate mesh with actual number of cross-sections (real + ghost)
    call mesh_initialise(msh, ncs)
    msh%has_ghost_cells = .true.
    
    ! Set actual cross-sections
    ics = 0
    do iseg = 1, nseg
      
        ! Upstream ghost cross-sections
        ics = ics + 1
        call crosssection_copy(cs(msh%seg(iseg)%first_cs), msh%cs(ics))
        msh%cs(ics)%ibase = 0
        ics = ics + 1
        call crosssection_copy(cs(msh%seg(iseg)%first_cs), msh%cs(ics))
        msh%cs(ics)%ibase = 0
        
        ! Interior cross-sections
        first_cs = ics + 1
        do i = msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
            ics = ics + 1
            call crosssection_copy(cs(i), msh%cs(ics))
        end do
        last_cs = ics
      
        ! Downstream ghost cross-sections
        ics = ics + 1
        call crosssection_copy(cs(msh%seg(iseg)%last_cs), msh%cs(ics))
        msh%cs(ics)%ibase = 0
        ics = ics + 1
        call crosssection_copy(cs(msh%seg(iseg)%last_cs), msh%cs(ics))
        msh%cs(ics)%ibase = 0
        
        ! Update indices of first and last cross-sections
        msh%seg(iseg)%first_cs = first_cs
        msh%seg(iseg)%last_cs = last_cs
        
    end do
    
    ! Compute polynomials and cumulated area and perimeters
    call update_geometries(msh)
    
    ! Compute curvilinear abscissae
    if (version < 1.3) then
        call curvilinear_abscissae_from_coords(msh)
    else
        call finalise_curvilinear_abscissae(msh)
    end if
    
    ! Compute slopes
    call bathy_slopes(msh)

    ! Check mesh
    call mesh_check(msh)
    
    return
  
100 write(buffer, '(I32)') iline
    write(err_msg, '(4A)') "error in file ", filename, ", at line ", trim(adjustl(buffer))
    close(10)
    call f90wrap_abort(err_msg)
101 write(buffer, '(I32)') iline
    write(err_msg, '(5A)') "error in file ", filename, ", at line ", trim(adjustl(buffer)), &
                           " (indices of upstream segments must be lower than index of current segment)"
    close(10)
    call f90wrap_abort(err_msg)
102 write(buffer, '(I32)') iline
    write(err_msg, '(5A)') "error in file ", filename, ", at line ", trim(adjustl(buffer)), &
                           " (index of downstream segment must be greater than index of current segment)"
    close(10)
    call f90wrap_abort(err_msg)
  
end subroutine
