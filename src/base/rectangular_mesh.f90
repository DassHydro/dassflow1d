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
subroutine rectangular_mesh(x, bathy, widths, msh)
    use m_mesh
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Array of curvilinear abscissae
    real(rp), dimension(:), intent(in) :: x
    !> Array of bathymetry elevations
    real(rp), dimension(:), intent(in) :: bathy
    !> Array of widths
    real(rp), dimension(:), intent(in) :: widths
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
   
    ! Check sizes
    if (size(x) /= size(bathy)) call abort_solver("bathy and x must be of same sizes")
    if (size(x) /= size(widths)) call abort_solver("widths and x must be of same sizes")
    
    ! Initialise mesh
    call mesh_initialise(msh, size(x), 1)
    msh%has_ghost_cells = .false.
    
    ! Create cross-sections
    do ics = 1, size(x)
    
        ! Initialise cross-section
        shape_model = "linear"
        call crosssection_initialise(msh%cs(ics), 1, shape_model)
        
        ! Set cross-section geometry
        msh%cs(ics)%coord%x = x(ics)
        msh%cs(ics)%coord%y = 0.0_rp
        msh%cs(ics)%x = x(ics)
        msh%cs(ics)%bathy = bathy(ics)
        msh%cs(ics)%level_heights(1) = bathy(ics) + 1000.0_rp
        msh%cs(ics)%level_widths(1) = widths(ics)

        ! Set overbanks levels
        msh%cs(ics)%ob_levels(1) = 1
        msh%cs(ics)%ob_levels(2) = 1
        
    end do
    
    ! Setup single segment
    msh%seg(1)%first_cs = 1
    msh%seg(1)%last_cs = size(x)
    allocate(msh%seg(1)%us_seg(1))
    nbc = 2
    msh%seg(1)%us_seg(1) = -1
    msh%seg(1)%us_bc = 1
    msh%seg(1)%ds_seg = -1
    msh%seg(1)%ds_bc = 2

    
    ! Compute polynomials and cumulated area and perimeters
    call update_geometries(msh)
    call finalise_curvilinear_abscissae(msh)
    
    ! Compute slopes
    call bathy_slopes(msh)
  
end subroutine
