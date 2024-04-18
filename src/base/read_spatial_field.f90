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
!> \file read_spatial_field.f90
!> \brief This file includes the routine for reading spatial fields from a file

!> Read a spatial field from a file
subroutine read_spatial_field(msh, filename)
    use m_common
    use m_mesh, only: Mesh
    use m_mesh, only: SpatialField
    use m_mesh, only: strickler_type_constant
    use m_mesh, only: strickler_type_powerlaw_h
    use m_mesh, only: strickler_type_einstein
    use m_mesh, only: strickler_type_average_conveyance
    use m_mesh, only: strickler_type_debord
    use m_mesh, only: set_pointwise_spatial_field
    use m_mesh, only: set_piecewise_uniform_spatial_field
    use m_mesh, only: set_piecewise_linear_spatial_field
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Mesh
    type(Mesh), intent(inout) :: msh
    !> Path to the file
    character(len=*), intent(in) :: filename
    ! TODO
    !=LOCAL VARIABLES==================================================================================================!
    ! Index
    integer(ip) :: i
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of field
    integer(ip) :: ifield
    ! Index of line in mesh file
    integer(ip) :: iline
    ! Index of segment
    integer(ip) :: iseg
    ! String for reading whole line in mesh file
    character(len=256) :: line
    ! Interpolation type
    character(len=lchar) :: interp
    ! Definition type (pointwise=all cross-sections, patch=by patche, defined by curvilinear abscissae)
    character(len=12) :: definition
    ! Number of cross-sections
    integer(ip) :: ncs
    ! Number of fields
    integer(ip) :: nfields
    ! Number of patches
    integer(ip) :: npatch
    ! Temporary value
    real(rp) :: tmp
    ! Mesh version
    real(rp) :: version
!     ! Mesh version
!     type(SpatialField), pointer :: field
    ! Temporary array for reading curvilinear abscissae of patches
    real(rp), dimension(:), allocatable :: x
    ! Temporary array for reading data of patches
    real(rp), dimension(:, :), allocatable :: y
    !> Variable
    character(len=12) :: var
    
    character(len=lchar) :: err_msg
   
    ! Open file
    open(10, file=filename, status="old")
    
    ! Read header line
    iline = 1
    read(10, '(A)', err=100) line
    if (line(1:27) == "DASSFLOW-1D SPATIAL FIELD V") then
        ! Read mesh version
        read(line(28:len_trim(line)), *, err=100) version
        
        ! Read header (comment) line
        iline = iline + 1
        read(10, '(A)', err=100) line
    else
        version = 1.0
    end if
    
    ! Read global parameters
    iline = iline + 1
    if (version == 1.0) then
        read(10, *, err=100) var, definition
    else
        goto 101
    end if
        
    ! Read definitions for each segment
    do iseg = 1, size(msh%seg)
    
        ! Skip header line for current segment (# SEGMENT XXX)
        iline = iline + 1
        read(10, '(A)', err=100) line
        
        ! Allocate spatial field(s)
        if (var == "strickler" .or. var == "K") then
        
            if (msh%strickler_type_code == strickler_type_constant) then
                nfields = 1
            else if (msh%strickler_type_code == strickler_type_powerlaw_h) then
                nfields = 2
            else if (msh%strickler_type_code == strickler_type_einstein) then
                nfields = 3
            else if (msh%strickler_type_code == strickler_type_average_conveyance) then
                nfields = 3
            else if (msh%strickler_type_code == strickler_type_debord) then
                nfields = 3
            end if
            if (.not. allocated(msh%seg(iseg)%strickler_fields)) then
                allocate(msh%seg(iseg)%strickler_fields(nfields))
            else
                goto 102
            end if
        else if (var == "bathy") then
            ! TODO
        else
            goto 103
        end if
        
        if (definition == "pointwise") then
            
            ncs = msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs + 1
            if (allocated(x)) deallocate(x)
            if (allocated(y)) deallocate(y)
            allocate(y(nfields, ncs))
            
            ! Read field values
            do i = 1, ncs
                iline = iline + 1
                read(10, *, err=100) y(:, i)
            end do
            
            ! Setup fields
            if (var == "strickler" .or. var == "K") then
                do ifield = 1, nfields
                    call set_pointwise_spatial_field(msh, iseg, msh%seg(iseg)%strickler_fields(ifield), y(ifield, :))
                end do
            else if (var == "bathy") then
                call set_pointwise_spatial_field(msh, iseg, msh%seg(iseg)%bathy_field, y(1, :))
            end if
            
        else
            
            iline = iline + 1
            read (10, *, err=100) npatch, interp
            if (allocated(x)) deallocate(x)
            allocate(x(npatch))
            if (allocated(y)) deallocate(y)
            allocate(y(nfields, npatch))
        
            ! Set x and y values for interpolation
            do i = 1, npatch
                iline = iline + 1
                read(10, *, err=100) x(i), y(:, i)
            end do
            
            ! Setup fields
            if (trim(interp) == "uniform") then
                if (var == "strickler" .or. var == "K") then
                    do ifield = 1, nfields
                        call set_piecewise_uniform_spatial_field(msh, iseg, msh%seg(iseg)%strickler_fields(ifield), &
                                                                 x, y(ifield, :))
                    end do
                else if (var == "bathy") then
                    call set_piecewise_uniform_spatial_field(msh, iseg, msh%seg(iseg)%bathy_field, &
                                                             x, y(1, :))
                end if
            else if (trim(interp) == "linear") then
                if (var == "strickler" .or. var == "K") then
                    do ifield = 1, nfields
                        call set_piecewise_linear_spatial_field(msh, iseg, msh%seg(iseg)%strickler_fields(ifield), &
                                                                x, y(ifield, :))
                    end do
                else if (var == "bathy") then
                    call set_piecewise_linear_spatial_field(msh, iseg, msh%seg(iseg)%bathy_field, &
                                                            x, y(1, :))
                end if
            else
                goto 104
            end if
            
        end if
        
    end do
    
    close(10)
    
    if (var == "strickler" .or. var == "K") then
        call apply_strickler_fields(msh)
    else if (var == "bathy") then
        call apply_bathy_field(msh)
    end if

    
    return
  
100 write(buffer, '(I32)') iline
    write(err_msg, '(4A)') "error in file ", filename, ", at line ", trim(adjustl(buffer))
    close(10)
    call f90wrap_abort(err_msg)
101 write(buffer, '(F32.2)') version
    write(err_msg, '(4A)') "wrong version in file ", filename, ": ", trim(adjustl(buffer))
    close(10)
    call f90wrap_abort(err_msg)
102 write(buffer, '(I32)') iseg
    write(err_msg, '(3A)') "strickler fields for segment ", trim(adjustl(buffer)), " already allocated"
    close(10)
    call f90wrap_abort(err_msg)
103 write(buffer, '(I32)') iline
    write(err_msg, '(5A)') "wrong variable name (", trim(var), ") in file ", filename, ", at line, ", &
                           trim(adjustl(buffer))
    close(10)
    call f90wrap_abort(err_msg)
104 write(buffer, '(I32)') iline
    write(err_msg, '(5A)') "wrong patch interpolation (", trim(interp), ") in file ", filename, ", at line, ", &
                           trim(adjustl(buffer))
    close(10)
    call f90wrap_abort(err_msg)
  
end subroutine
