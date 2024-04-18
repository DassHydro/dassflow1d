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
!> @file m_sw_mono.f90
!> @brief This file includes the m_model module.

!> @brief Module m_sw_mono.
module m_sw_mono
    use m_common
    use m_linear_algebra
    use m_mesh, only: Mesh
    use m_obs
!     use m_time_screen                                                                                              !NOADJ
    implicit none
#ifndef CPP_ADJ
    include "mpif.h"
    include 'dmumps_struc.h'
#endif

   !===================================================================================================================!
   !  Discrete Model Unknowns Structure
   !===================================================================================================================!

    !> @brief Structure of unknowns.
    !> @details This structure includes all the unknows of models (A, Q, h, qlat) and the terms for compuing the missfit
    type Unknowns

        !> Flow areas
        real(rp), dimension(:), allocatable :: A
        !> Discharges
        real(rp), dimension(:), allocatable :: Q
        !> Depths
        real(rp), dimension(:), allocatable :: h
        !> Lateral inflows
        real(rp), dimension(:,:), allocatable :: qlat
        !> Friction term values
        real(rp), dimension(:), allocatable :: Sf
        !> Gravity term values
        real(rp), dimension(:), allocatable :: Sg

    end type

    !> Structure of implicit matrix
    type ImplicitMatrix

        !> Sizes
        integer(ip), dimension(2) :: sizes
        !> Segments offsets
        integer(ip), dimension(:), allocatable :: seg_offsets
        !> GA values
        real(rp), dimension(:), allocatable :: GA
        !> GB values
        real(rp), dimension(:), allocatable :: GB
        !> GC values
        real(rp), dimension(:), allocatable :: GC
        !> GD values
        real(rp), dimension(:), allocatable :: GD
        !> GE values
        real(rp), dimension(:), allocatable :: GE
        !> GF values
        real(rp), dimension(:), allocatable :: GF
        !> CR values
        real(rp), dimension(:), allocatable :: CR
        !> CS values
        real(rp), dimension(:), allocatable :: CS
        !> CT values
        real(rp), dimension(:), allocatable :: CT
        
        ! TA values (condensed relations x1 = TA.xn + TB)
        real(rp), dimension(:), allocatable :: TA1
        real(rp), dimension(:), allocatable :: TA2
        real(rp), dimension(:), allocatable :: TA3
        real(rp), dimension(:), allocatable :: TA4
        ! TB values (condensed relations x1 = TA.xn + TB)
        real(rp), dimension(:), allocatable :: TB1
        real(rp), dimension(:), allocatable :: TB2

#ifndef CPP_ADJ     
        ! MUMPS System (condensed system)
        type(DMUMPS_STRUC) :: sys
#endif
        real(rp), dimension(:), allocatable :: ANZ
        real(rp), dimension(:), allocatable :: RHS

    end type

    !> Structure of timeseries.
    type Timeseries
        real(rp), dimension(:), allocatable  ::  t  !< Time array.
        real(rp), dimension(:), allocatable  ::  y  !< Flow array.
    end type

    !> Structure of boundary condition
    type BoundaryCondition

        !> ID
        character(len=16) :: id
        !> Timeseries
        type(Timeseries) :: ts

    end type

    !> Structure of inflow condition
    type InflowCondition

        !> Index of segment
        integer(ip) :: iseg
        !> Index of edge (between cross-sections [ie, ie+1])
        integer(ip) :: ie
        !> Length of edge (q=ts%y * ledge)
        real(rp) :: ledge
        !> Timeseries
        type(Timeseries) :: ts

    end type
    
    type Results
    
        real(rp), dimension(:), allocatable :: t
        real(rp), dimension(:, :), allocatable :: q
        real(rp), dimension(:, :), allocatable :: h
        real(rp), dimension(:, :), allocatable :: a
    
    end type
    
    type Probe
        integer(ip) :: ics
        integer(ip) :: it
        real(rp) :: dt
        real(rp), dimension(:), allocatable :: t
        real(rp), dimension(:, :), allocatable :: data
        character*1, dimension(:), allocatable :: variables
    end type Probe
        
    type Model

        !> Flag to enable multiple steady states rather than unsteady
        logical :: steady_states
        !> Flag to enable discharge in estimations
        logical :: discharge_estimation
        !> True if initialisation is done externally
        logical :: external_initialisation
        !> True to write observations status
        logical :: observations_status
        !> True to write results at resampled cross-sections in outputs
        logical :: output_resampled_cs
        !> True to print progress on screen
        logical :: print_progress
        !> True to disable standard output
        logical :: disable_stdout
        !> Status code
        !> Possible values are:
        !> - 10 : NaN value(s) detected at downstream boundary during standard step.
        !> - 11 : NaN value(s) detected inside a segment during standard step.
        !> - 20 : NaN value(s) detected on hydraulic variables (after geometry computation).
        !> - 21 : NaN value(s) detected on h or Q.
        !> - 22 : Negative Debord operand.
        integer(ip) :: status
        !> Warning counters (1=number of heps correction encountered, 2=TBD, 3=TBD, 4=TBD, 5=TBD)
        integer(ip), dimension(5) :: warning_counters
        !> Internal counters (1=overbanks flow, 2=left overbank flow, 3=right overbank flow, 4=TBD, 5=TBD)
        integer(ip), dimension(5) :: internal_counters
        !> Current time
        real(rp) :: tc
        !> End time
        real(rp) :: te
        !> Start time
        real(rp) :: ts
        !> Timestep
        real(rp) :: dt
        !> Timestep for printing informations on standard output
        real(rp) :: dtout
        !> Cutoff for gradients of wetted perimeter w.r.t. height
        real(rp) :: dpdz_cutoff
        !> Froude threshold for LPI method
        real(rp) :: frLPI
        !> Gravity acceleration (m2/s)
        real(rp) :: gravity
        !> Depth threshold
        real(rp) :: heps
        !> Discharge threshold
        real(rp) :: qeps
        !> Parameter for LPI method
        real(rp) :: mLPI
        !> Implicit coefficient for Preissmann scheme
        real(rp) :: theta_preissmann
        !> Regularization coefficient
!         real(rp) :: alpha_reg
        real(rp) :: gamma_reg
        !> Observation cost
        real(rp) :: cost_obs
        !> Regularization cost
        real(rp) :: cost_reg
        !> Output file
        character(len=128) :: output_file
        !> Computational scheme
        character(len=32) :: scheme
        !>scale for the implicit diffusive wave
        character(len=16) :: scale_idw
        !> Bounday conditions
        type(BoundaryCondition), dimension(:), allocatable :: bc
        !> Number of inflow conditions
        integer(ip) :: nic
        !> Inflow conditions
        type(InflowCondition), dimension(:), allocatable :: ic
        !> Mesh
        type(Mesh) :: msh
        !> large grid
        type(Mesh), pointer :: large_grid
        !> Unknowns
        type(Unknowns) :: dof
        !> Implicit arrays
        type(ImplicitMatrix) :: imp
        !> Results
        type(Results) :: res
        !> Probes
        type(Probe), dimension(:), allocatable :: probes
      
    end type


#ifndef CPP_ADJ
    interface
        subroutine generate_observations(mdl, obs)
            import Model
            import Observations
            implicit none
            type(Model), intent(inout) :: mdl
            type(Observations), intent(inout), optional :: obs
        end subroutine
        subroutine standard_step(mdl, eps, itermax)
            import ip
            import Model
            import rp
            implicit none
            type(Model), intent(inout) :: mdl
            real(rp), intent(in) :: eps
            integer(ip), intent(in) :: itermax
        end subroutine
        subroutine preissmann_timestep(mdl, msh, imp, dof, status)
            import ip
            import ImplicitMatrix
            import Mesh
            import Model
            import Unknowns
            implicit none
            type(Model), intent(in) :: mdl
            type(Mesh), intent(inout) :: msh
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
            integer(ip) :: status
        end subroutine
	subroutine implicit_diffusive_wave(mdl, msh, large_grid, imp, dof, status)
            import ip
            import ImplicitMatrix
            import Mesh
            import Model
            import Unknowns
            implicit none
            type(Model), intent(in) :: mdl
            type(Mesh), intent(inout) :: msh
            type(Mesh), intent(inout) :: large_grid
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
            integer(ip) :: status
        end subroutine
!         subroutine steady_state(mdl, obs)
!             import Model
!             import Observations
!             implicit none
!             type(Model), intent(inout) :: mdl
!             type(Observations), intent(inout), optional :: obs
!         end subroutine
!         subroutine steady_states_loop(mdl, obs)
!             import Model
!             import Observations
!             implicit none
!             type(Model), intent(inout) :: mdl
!             type(Observations), intent(inout), optional :: obs
!         end subroutine
        subroutine time_loop(mdl, obs)
            import Model
            import Observations
            implicit none
            type(Model), intent(inout) :: mdl
            type(Observations), intent(inout), optional :: obs
        end subroutine
!         subroutine free_surface_slopes_segment(mdl, iseg, h, slopes)
!             import ip
!             import rp
!             import Model
!             import Unknowns
!             implicit none
!             type(Model), intent(in) :: mdl
!             integer(ip), intent(in) :: iseg
!             real(rp), dimension(:), intent(in) :: h
!             real(rp), dimension(:), intent(inout) :: slopes
!         end subroutine
    end interface
#endif

      
#ifndef CPP_ADJ
    contains
    
    !==================================================================================================================!
    !  Methods for the Unknows type
    !==================================================================================================================!

     !> @brief Initialise the unknowns of the model.
    subroutine unknowns_initialise(dof, msh)
        implicit none
        !=ARGUMENTS====================================================================================================!
        !> Unknowns
        type(Unknowns), intent(inout) :: dof
        !> Mesh
        type(Mesh), intent(in ) :: msh

        ! Allocate arrays
        allocate(dof%Q(msh%ncs))
        allocate(dof%A(msh%ncs))
        allocate(dof%h(msh%ncs))
        allocate(dof%sg(msh%ncs))
        allocate(dof%sf(msh%ncs))


        ! Set initial values
        dof%Q(:)  =  0._rp
        dof%A(:)  =  0._rp
        dof%h(:)  =  0._rp
        dof%Sg(:)  =  0._rp
        dof%Sf(:) =  0._rp
        
    end subroutine


    subroutine unknowns_finalise(dof)
        implicit none
        !=ARGUMENTS====================================================================================================!
        !> Unknowns
        type(Unknowns), intent(inout) :: dof
        
        if (allocated(dof%Q)) deallocate(dof%Q)
        if (allocated(dof%A)) deallocate(dof%A)
        if (allocated(dof%h)) deallocate(dof%h)
        if (allocated(dof%Sf)) deallocate(dof%Sf)
        if (allocated(dof%Sg)) deallocate(dof%Sg)
        
    end subroutine    
    
    !==================================================================================================================!
    !  Methods for the ImplicitMatrix type
    !==================================================================================================================!
    
    subroutine implicitmatrix_initialise(imp, mdl, msh)
        implicit none
        type(Model), intent(in) :: mdl
        type(Mesh), intent(in) :: msh
        type(ImplicitMatrix), intent(inout) :: imp
        
        ! Index of segment
        integer(ip) :: iseg
        ! Number of rows
        integer(ip) :: nrows
        ! Number of boundary conditions
        integer(ip) :: nbc
        
        ! Allocate offsets array
        allocate(imp%seg_offsets(msh%nseg))
        
        ! Compute offsets and total number of rows
        nrows = 0
        do iseg = 1, msh%nseg
        
            imp%seg_offsets(iseg) = nrows
            nrows = nrows + msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs + 1

        end do
        
        ! Allocate coefficients arrays
        allocate(imp%GA(nrows))
        allocate(imp%GB(nrows))
        allocate(imp%GC(nrows))
        allocate(imp%GD(nrows))
        allocate(imp%GE(nrows))
        allocate(imp%GF(nrows))
        allocate(imp%CR(nrows))
        allocate(imp%CS(nrows))
        allocate(imp%CT(nrows))
        
        ! Allocate coefficients arrays for condensed relations
        allocate(imp%TA1(msh%nseg))
        allocate(imp%TA2(msh%nseg))
        allocate(imp%TA3(msh%nseg))
        allocate(imp%TA4(msh%nseg))
        allocate(imp%TB1(msh%nseg))
        allocate(imp%TB2(msh%nseg))
        
        call init_implicit_reduced_system(mdl, msh, imp)
        
    end subroutine
    
    
    subroutine implicitmatrix_finalise(imp)
        implicit none
        type(ImplicitMatrix), intent(inout) :: imp
        
        ! Deallocate arrays
        if (allocated(imp%seg_offsets)) deallocate(imp%seg_offsets)
        if (allocated(imp%GA)) deallocate(imp%GA)
        if (allocated(imp%GB)) deallocate(imp%GB)
        if (allocated(imp%GC)) deallocate(imp%GC)
        if (allocated(imp%GD)) deallocate(imp%GD)
        if (allocated(imp%GE)) deallocate(imp%GE)
        if (allocated(imp%GF)) deallocate(imp%GF)
        if (allocated(imp%CR)) deallocate(imp%CR)
        if (allocated(imp%CS)) deallocate(imp%CS)
        if (allocated(imp%CT)) deallocate(imp%CT)
        if (allocated(imp%TA1)) deallocate(imp%TA1)
        if (allocated(imp%TA2)) deallocate(imp%TA2)
        if (allocated(imp%TA3)) deallocate(imp%TA3)
        if (allocated(imp%TA4)) deallocate(imp%TA4)
        if (allocated(imp%TB1)) deallocate(imp%TB1)
        if (allocated(imp%TB2)) deallocate(imp%TB2)
        
    end subroutine
    
    !==================================================================================================================!
    !  Methods for the Timeseries type
    !==================================================================================================================!
    
    subroutine timeseries_initialise(ts, nt, y0)
        implicit none
        type(Timeseries), intent(out) :: ts
        integer(ip), intent(in) :: nt
        real(rp), optional :: y0
        
        if (nt == 0) then
            allocate(ts%y(1))
        else
            allocate(ts%t(nt))
            allocate(ts%y(nt))
        end if
        if (present(y0)) then
            ts%y(:) = y0
        else
            ts%y(:) = 0.0_rp
        end if
        
    end subroutine
    
    
    subroutine timeseries_finalise(ts)
        implicit none
        type(Timeseries), intent(inout) :: ts
        
        ! Deallocate array
        if (allocated(ts%t)) deallocate(ts%t)
        if (allocated(ts%y)) deallocate(ts%y)
        
    end subroutine
    
    
    subroutine reset_timeseries(ts, t, y)
        implicit none
        type(Timeseries), intent(inout) :: ts
        real(rp), dimension(:), intent(in) :: t
        real(rp), dimension(:), intent(in) :: y
        
        ! CHECK-UP
        if(size(y) /= size(t)) then
            call abort_solver("t and y must be of same size")
        end if
        
        ! Reallocate arrays
        if (allocated(ts%t)) deallocate(ts%t)
        allocate(ts%t(size(t)))
        if (allocated(ts%y)) deallocate(ts%y)
        allocate(ts%y(size(t)))
        
        ! Set arrays values
        ts%t(:) = t(:)
        ts%y(:) = y(:)
        
    end subroutine
    
    !==================================================================================================================!
    !  Methods for the BoundaryCondition type
    !==================================================================================================================!
    
    subroutine set_timeseries(bc, t, y)
        implicit none
        type(BoundaryCondition), intent(inout) :: bc
        real(rp), dimension(:), intent(in) :: t
        real(rp), dimension(:), intent(in) :: y
        
        integer(ip) :: i

        if (size(t) /= size(y)) then
            call f90wrap_abort("t and y must be of same size")
        end if
        
        ! CHECK-UP
        do i = 2, size(t)
            if (t(i) <= t(i-1)) call f90wrap_abort("t values must be in ascending order")
        end do
        
        if (allocated(bc%ts%t)) deallocate(bc%ts%t)
        if (allocated(bc%ts%y)) deallocate(bc%ts%y)
        
        allocate(bc%ts%t(size(t)))
        allocate(bc%ts%y(size(t)))
        bc%ts%t(:) = t
        bc%ts%y(:) = y
        
    end subroutine
    
    
    subroutine set_rating_curve(bc, z, q)
        implicit none
        type(BoundaryCondition), intent(inout) :: bc
        real(rp), dimension(:), intent(in) :: z
        real(rp), dimension(:), intent(in) :: q
        
        integer(ip) :: i

        if (size(z) /= size(q)) then
            call f90wrap_abort("z and q must be of same size")
        end if
        
        ! CHECK-UP
        do i = 2, size(z)
            if (z(i) <= z(i-1)) call f90wrap_abort("z values must be in ascending order")
            if (q(i) <= q(i-1)) call f90wrap_abort("q values must be in ascending order")
        end do
        
        if (allocated(bc%ts%t)) deallocate(bc%ts%t)
        if (allocated(bc%ts%y)) deallocate(bc%ts%y)
        
        allocate(bc%ts%t(size(z)))
        allocate(bc%ts%y(size(z)))
        bc%ts%t(:) = z
        bc%ts%y(:) = q
        
    end subroutine
    
    !==================================================================================================================!
    !  Methods for the InflowCondition type
    !==================================================================================================================!
    
    subroutine inflowcondition_finalise(ic)
        implicit none
        type(InflowCondition), intent(inout) :: ic

        call timeseries_finalise(ic%ts)
        
    end subroutine
    
    subroutine copy_inflow_condition(src, dst)
        implicit none
        type(InflowCondition), intent(in) :: src
        type(InflowCondition), intent(inout) :: dst

        dst%iseg = src%iseg
        dst%ie = src%ie
        dst%ledge= src%ledge
        if (allocated(src%ts%t)) then
            allocate(dst%ts%t(size(src%ts%t)))
            dst%ts%t(:) = src%ts%t(:)
            allocate(dst%ts%y(size(src%ts%y)))
            dst%ts%y(:) = src%ts%y(:)
        end if
        
    end subroutine
    
    subroutine move_inflow_condition(src, dst)
        implicit none
        type(InflowCondition), intent(inout) :: src
        type(InflowCondition), intent(inout) :: dst

        call copy_inflow_condition(src, dst)
        call inflowcondition_finalise(src)
        
    end subroutine
    
    function get_discharge_timeseries(ic) result(ts)
        implicit none
        type(InflowCondition), intent(in) :: ic
        type(TimeSeries) :: ts

        if (.not. allocated(ic%ts%y)) call f90wrap_abort("Inflow timeseries not set")

        allocate(ts%t(size(ic%ts%t)))
        allocate(ts%y(size(ic%ts%t)))
        ts%t(:) = ic%ts%t(:)
        ts%y(:) = ic%ts%y(:) * ic%ledge
        
    end function

    
    subroutine reset_discharge_timeseries(ic, t, q)
        implicit none
        type(InflowCondition), intent(inout) :: ic
        real(rp), dimension(:) :: t
        real(rp), dimension(:) :: q

        ! CHECK-UP
        if (size(t) /= size(q)) then
            call f90wrap_abort("t and q must be of same size")
        end if

        ! Reset timeseries
        call reset_timeseries(ic%ts, t, q / ic%ledge)
        
    end subroutine
    
    !==================================================================================================================!
    !  Methods for the Probe type
    !==================================================================================================================!
    
    subroutine probe_initialise(prb, msh, ics, dt, variables, base_cs)
        implicit none
        type(Probe), intent(out) :: prb
        type(Mesh), intent(in), target :: msh
        integer(ip), intent(in) :: ics
        real(rp), intent(in) :: dt
        character(len=*), intent(in) :: variables
        logical, intent(in), optional :: base_cs
        
        integer(ip) :: i1
        integer(ip) :: ivar
        logical :: use_base_cs

        if (present(base_cs)) then
            use_base_cs = base_cs
        else
            use_base_cs = .false.
        end if

        ! Setup timestep
        prb%dt = dt
        
        ! Setup array of variables identifiers
        allocate(prb%variables(len(variables)))
        do ivar = 1, len(variables)
            prb%variables(ivar) = variables(ivar:ivar)
        end do
        
        ! Retrieve index of base cross-section
        if (use_base_cs) then
            do i1 = 1, msh%ncs
                if (msh%cs(i1)%ibase == ics) then
                    prb%ics = i1
                end if
            end do
            ! TODO Check base ics found
        else
            ! TODO Check ics does not overflow
            prb%ics = ics
        end if
        
    end subroutine

    
    subroutine probe_finalise(prb)
        implicit none
        type(Probe), intent(inout) :: prb
        
        if (allocated(prb%variables)) deallocate(prb%variables)
        if (allocated(prb%t)) deallocate(prb%t)
        if (allocated(prb%data)) deallocate(prb%data)
        
    end subroutine

    
    subroutine copy_probe(src, dst)
        implicit none
        type(Probe), intent(in) :: src
        type(Probe), intent(inout) :: dst
        
        ! Clear destination probe
        call probe_finalise(dst)
        
        ! Copy cross-section index, timestep and variables
        dst%ics = src%ics
        dst%dt = src%dt
        if (allocated(src%variables)) then
            allocate(dst%variables(size(src%variables)))
            dst%variables(:) = src%variables(:)
        end if
        
    end subroutine
    
    
    subroutine model_initialise(mdl, msh, scheme)
        implicit none
        type(Mesh), intent(in), target :: msh
        type(Model), intent(out) :: mdl
        character(len=*), intent(in), optional :: scheme
        
        ! Index of boundary condition
        integer(ip) :: ibc
        ! Index of segment
        integer(ip) :: iseg
        ! Number of boundary conditions
        integer(ip) :: nbc

        ! Copy mesh
        call mesh_copy(msh, mdl%msh)

        ! Set large grid default pointer
        mdl%large_grid => NULL()       
        
        ! Compute number of boundary conditions
        nbc = 0
        do iseg = 1, msh%nseg
        
            if (msh%seg(iseg)%us_bc > 0) then
                nbc = max(nbc, msh%seg(iseg)%us_bc)
            end if
            if (msh%seg(iseg)%ds_bc > 0) then
                nbc = max(nbc, msh%seg(iseg)%ds_bc)
            end if

        end do
        
        ! Allocate and initialise boundary conditions
        ! TODO check for already allocated bcs
        allocate(mdl%bc(nbc))
        do ibc = 1, nbc
            mdl%bc(ibc)%id = "undefined"
        end do
        
        ! Set default values
        mdl%steady_states = .false.
        mdl%discharge_estimation = .false.
        mdl%external_initialisation = .false.
        mdl%print_progress = .true.
        mdl%output_resampled_cs = .false.
        mdl%observations_status = .false.
        mdl%disable_stdout = .false.
        mdl%output_file = ""
        mdl%tc = 0.0
        mdl%ts = 0.0
        mdl%te = 1.0
        mdl%dt = 1.0
        mdl%dtout = 1.0
        mdl%heps = 0.01
        mdl%qeps = 0.001
        mdl%dpdz_cutoff = huge(1.0_rp)
        mdl%frLPI = 0.7
        mdl%mLPI = 10.0
        mdl%gravity = 9.81
        mdl%theta_preissmann = 0.7
!         mdl%alpha_reg = 0.0
        mdl%gamma_reg = -1.0
        mdl%scale_idw = "low"
        
        ! Allocate unknowns
        call unknowns_initialise(mdl%dof, mdl%msh)
        
        mdl%scheme = "undefined"
        mdl%nic = 0
        
    end subroutine
    
    subroutine set_large_grid(mdl,large_grid) 

	implicit none
        type(Model), intent(inout) :: mdl
        type(Mesh), intent(in), target :: large_grid
        
        mdl%large_grid => large_grid

    end subroutine    

    subroutine add_inflow_condition(mdl, iseg, coords, t, q)
        implicit none
        type(Model) :: mdl
        integer(ip) :: iseg
        real(rp), dimension(2) :: coords
        real(rp), dimension(:) :: t
        real(rp), dimension(:) :: q
        
        ! Iterator
        integer(ip) :: i
        ! Index of cross-section
        integer(ip) :: ics
        ! Index of edge
        integer(ip) :: ie
        ! Index of inflow condition
        integer(ip) :: iic
        ! Number of inflow conditions before adding new inflow condition
        integer(ip) :: size_old
        ! Squares of distances
        real(rp), dimension(2) :: dist2
        ! Minimum square of distance
        real(rp) :: dist2min
        ! Temporary arry for copying existing inflow conditions
        type(InflowCondition), dimension(:), allocatable :: old_ic
        
        ! CHECK-UP
        if (iseg < 0 .or. iseg > size(mdl%msh%seg)) then
            call f90wrap_abort("index of segment out of bounds")
        end if
        if (size(t) /= size(q)) then
            call f90wrap_abort("t and q must be of same size")
        end if

        ! Locate closest cross-section
        ie = mdl%msh%seg(iseg+1)%first_cs
        ics = ie
        dist2min = (coords(1) - mdl%msh%cs(ics)%coord%x)**2 + (coords(2) - mdl%msh%cs(ics)%coord%y)**2
        do ics = mdl%msh%seg(iseg+1)%first_cs+1, mdl%msh%seg(iseg+1)%last_cs
            dist2(1) = (coords(1) - mdl%msh%cs(ics)%coord%x)**2 + (coords(2) - mdl%msh%cs(ics)%coord%y)**2
            if (dist2min > dist2(1)) then
                dist2min = dist2(1)
                ie = ics
            end if
        end do
        
        ! Update index of edge
        if (ie > mdl%msh%seg(iseg+1)%first_cs .and. ie < mdl%msh%seg(iseg+1)%last_cs) then
            dist2(1) = (coords(1) - mdl%msh%cs(ie-1)%coord%x)**2 + (coords(2) - mdl%msh%cs(ie-1)%coord%y)**2
            dist2(2) = (coords(1) - mdl%msh%cs(ie+1)%coord%x)**2 + (coords(2) - mdl%msh%cs(ie+1)%coord%y)**2
            if (dist2(1) < dist2(2)) ie = ie - 1
        end if
        
        ! Check that inflow conditions are set in order (increasing iseg, ics)
        if (allocated(mdl%ic) .and. mdl%nic > 0) then
            size_old = size(mdl%ic)
            if (iseg+1 < mdl%ic(size_old)%iseg) then
                call f90wrap_abort("inflow conditions must be added with increasing indices of segments")
            end if
            if (iseg+1 == mdl%ic(size_old)%iseg .and. ie < mdl%ic(size_old)%ie) then
                call f90wrap_abort("inflow conditions must be added with increasing indices of edges in a segment")
            end if
        end if
        
        
        ! Resize array of inflow conditions
        if (allocated(mdl%ic)) then
            size_old = size(mdl%ic)
            allocate(old_ic(size(mdl%ic)))
            do i = 1, size(mdl%ic)
                call move_inflow_condition(mdl%ic(i), old_ic(i))
!                 old_ic(i)%iseg = mdl%ic(i)%iseg
!                 old_ic(i)%ie = mdl%ic(i)%ie
!                 old_ic(i)%ledge= mdl%ic(i)%ledge
!                 allocate(old_ic(i)%ts%t(size(mdl%ic(i)%ts%t)))
!                 old_ic(i)%ts%t(:) = mdl%ic(i)%ts%t(:)
!                 allocate(old_ic(i)%ts%y(size(mdl%ic(i)%ts%y)))
!                 old_ic(i)%ts%y(:) = mdl%ic(i)%ts%y(:)
!                 deallocate(mdl%ic(i)%ts%t)
!                 deallocate(mdl%ic(i)%ts%y)
            end do
            deallocate(mdl%ic)
        else
            size_old = 0
        end if
        allocate(mdl%ic(size_old+1))
        if (size_old > 0) then
            do i = 1, size_old
                call move_inflow_condition(old_ic(i), mdl%ic(i))
!                 mdl%ic(i)%iseg = old_ic(i)%iseg
!                 mdl%ic(i)%ie = old_ic(i)%ie
!                 mdl%ic(i)%ledge = old_ic(i)%ledge
!                 allocate(mdl%ic(i)%ts%t(size(old_ic(i)%ts%t)))
!                 mdl%ic(i)%ts%t(:) = old_ic(i)%ts%t(:)
!                 allocate(mdl%ic(i)%ts%y(size(old_ic(i)%ts%t)))
!                 mdl%ic(i)%ts%y(:) = old_ic(i)%ts%y(:)
!                 deallocate(old_ic(i)%ts%t)
!                 deallocate(old_ic(i)%ts%y)
            end do
            deallocate(old_ic)
        end if
        
        ! Set new inflow condition
        mdl%ic(size_old+1)%iseg = iseg + 1
        mdl%ic(size_old+1)%ie = ie
        mdl%ic(size_old+1)%ledge = mdl%msh%cs(ie)%deltademi
        allocate(mdl%ic(size_old+1)%ts%t(size(t)))
        mdl%ic(size_old+1)%ts%t(:) = t(:)
        allocate(mdl%ic(size_old+1)%ts%y(size(q)))
        mdl%ic(size_old+1)%ts%y(:) = q(:)
        
        mdl%nic = size(mdl%ic)
        
    end subroutine
    
    
    subroutine set_scheme(mdl, scheme)
        implicit none
        type(Model), intent(inout) :: mdl
        character(len=*), intent(in), optional :: scheme
        
        ! Index of boundary condition
        integer(ip) :: ibc
        ! Index of segment
        integer(ip) :: iseg
        ! Number of boundary conditions
        integer(ip) :: nbc
    
        ! CHECK-UP
        do ibc = 1, size(mdl%bc)
        
            if (mdl%bc(ibc)%id == "undefined") then
                call abort_solver("All boundary conditions must be set before calling 'set_scheme'")
            end if
        
        end do
        ! Allocate implicit arrays
        if (present(scheme)) then
            if (scheme == "preissmann") then
                mdl%scheme = "preissmann"
                call implicitmatrix_initialise(mdl%imp, mdl, mdl%msh)
                
            end if
	    if (scheme == "implicit_diffusive_wave") then
                mdl%scheme = "implicit_diffusive_wave"
                call implicitmatrix_initialise(mdl%imp, mdl, mdl%msh)
            end if
            
        else
            mdl%scheme = "implicit_diffusive_wave"
            !mdl%scheme = "preissmann"
            call implicitmatrix_initialise(mdl%imp, mdl, mdl%msh)
        
        end if
        
    end subroutine
    
    subroutine add_probe(mdl, ics, dt, variables, base_cs)
        implicit none
        type(Model), intent(inout) :: mdl
        integer(ip), intent(in) :: ics
        real(rp), intent(in) :: dt
        character(len=*), intent(in) :: variables
        logical, intent(in), optional :: base_cs
        
        integer(ip) :: i1
        integer(ip) :: iprobe
        integer(ip) :: ivar
        integer(ip) :: nprobes
        logical :: use_base_cs
        type(Probe), dimension(:), allocatable :: probes_tmp
!         
!         do ivar = 1, len(variables)
!             print *, ivar, variables(ivar:ivar)
!         end do

        if (present(base_cs)) then
            use_base_cs = base_cs
        else
            use_base_cs = .false.
        end if
        
        ! Store array of already defined probes
        if (allocated(mdl%probes)) then
        
            allocate(probes_tmp(size(mdl%probes)))
            do iprobe = 1, size(mdl%probes)
                call copy_probe(mdl%probes(iprobe), probes_tmp(iprobe))
                call probe_finalise(mdl%probes(iprobe))
            end do
            nprobes = size(mdl%probes)
            
        else
        
            nprobes = 0
        
        end if
        
        ! Reallocate array of probes
        if (allocated(mdl%probes)) deallocate(mdl%probes)
        allocate(mdl%probes(nprobes+1))

        ! Copy back already defined probes
        if (nprobes > 0) then
        
            do iprobe = 1, nprobes
                call copy_probe(probes_tmp(iprobe), mdl%probes(iprobe))
                call probe_finalise(probes_tmp(iprobe))
            end do
        
        end if
        
        ! Setup new probe
        call probe_initialise(mdl%probes(nprobes+1), mdl%msh, ics, dt, variables, use_base_cs)
        
    end subroutine
    
    
    subroutine add_probe_using_coords(mdl, coords, dt, variables, base_cs)
        implicit none
        type(Model), intent(inout) :: mdl
        real(rp), dimension(2), intent(in) :: coords
        real(rp), intent(in) :: dt
        character(len=*), intent(in) :: variables
        logical, intent(in), optional :: base_cs
        
        ! Iterators
        integer(ip) :: i1
        integer(ip) :: mindistindex
        real(rp) :: dist
        real(rp) :: mindist
        integer(ip) :: iprobe
        integer(ip) :: ivar
        integer(ip) :: nprobes
        logical :: use_base_cs
        type(Probe), dimension(:), allocatable :: probes_tmp
        
        if (present(base_cs)) then
            use_base_cs = base_cs
        else
            use_base_cs = .false.
        end if

        ! Store array of already defined probes
        if (allocated(mdl%probes)) then
        
            allocate(probes_tmp(size(mdl%probes)))
            do iprobe = 1, size(mdl%probes)
                call copy_probe(mdl%probes(iprobe), probes_tmp(iprobe))
                call probe_finalise(mdl%probes(iprobe))
            end do
            nprobes = size(mdl%probes)
            
        else
        
            nprobes = 0
        
        end if
        
        ! Reallocate array of probes
        if (allocated(mdl%probes)) deallocate(mdl%probes)
        allocate(mdl%probes(nprobes+1))

        ! Copy back already defined probes
        if (nprobes > 0) then
        
            do iprobe = 1, nprobes
                call copy_probe(probes_tmp(iprobe), mdl%probes(iprobe))
                call probe_finalise(probes_tmp(iprobe))
            end do
        
        end if
        
        ! Find index of cross-section
        do i1 = 1, mdl%msh%ncs
            dist = sqrt((coords(1) - mdl%msh%cs(i1)%coord%x)**2 + &
                        (coords(2) - mdl%msh%cs(i1)%coord%y)**2)
            if (i1 == 1) then
                mindist = dist
                mindistindex = i1
            else if (dist < mindist) then
                mindist = dist
                mindistindex = i1
            end if
        end do
        
        ! Setup new probe
        call probe_initialise(mdl%probes(nprobes+1), mdl%msh, mindistindex, dt, variables, use_base_cs)
        
    end subroutine add_probe_using_coords

    
    subroutine reset_probes(mdl)
        implicit none
        type(Model), intent(inout) :: mdl
        
        integer(ip) :: iprobe
        integer(ip) :: it
        integer(ip) :: nt
        
        if (.not. allocated(mdl%probes)) return
        
        do iprobe = 1, size(mdl%probes)
        
            nt = floor((mdl%te - mdl%ts) / mdl%probes(iprobe)%dt)
            
            ! (Re)allocate time and data arrays
            if (allocated(mdl%probes(iprobe)%t)) deallocate(mdl%probes(iprobe)%t)
            allocate(mdl%probes(iprobe)%t(nt))
            if (allocated(mdl%probes(iprobe)%data)) deallocate(mdl%probes(iprobe)%data)
            allocate(mdl%probes(iprobe)%data(size(mdl%probes(iprobe)%variables), nt))
            
            ! Set time values
            do it = 1, nt
                mdl%probes(iprobe)%t(it) = mdl%ts + mdl%probes(iprobe)%dt * it
            end do
            
            ! Set NA data values
            mdl%probes(iprobe)%data(:, :) = -9999.0
            
            ! Reset counter
            mdl%probes(iprobe)%it = 1
        
        end do
    
    end subroutine

    
    subroutine run_steady(mdl, obs)
        implicit none
        type(Model), intent(inout) :: mdl
        type(Observations), intent(inout), optional :: obs
        
        call steady_state(mdl, obs)
        
    end subroutine

    
    subroutine run_steady_batch(mdl, obs)
        implicit none
        type(Model), intent(inout) :: mdl
        type(Observations), intent(inout), optional :: obs
        
        call steady_states_loop(mdl, obs)
        
    end subroutine

    
    subroutine run_unsteady(mdl, obs)
        implicit none
        type(Model), intent(inout) :: mdl
        type(Observations), intent(inout), optional :: obs
        
        call time_loop(mdl, obs)
        
    end subroutine
    
    
    subroutine model_finalise(mdl)
        implicit none
        type(Model), intent(inout) :: mdl
        
        integer(ip) :: i
        
        ! Deallocate boundary conditions
        if (allocated(mdl%bc)) then
            do i = 1, size(mdl%bc)
                call timeseries_finalise(mdl%bc(i)%ts)
            end do
            deallocate(mdl%bc)
        end if
        
        ! Deallocate inflow conditions
        if (allocated(mdl%ic)) then
            do i = 1, size(mdl%ic)
                call timeseries_finalise(mdl%ic(i)%ts)
            end do
            deallocate(mdl%ic)
        end if
        
        ! Deallocate probes
        if (allocated(mdl%probes)) then
            do i = 1, size(mdl%probes)
                call probe_finalise(mdl%probes(i))
            end do
            deallocate(mdl%probes)
        end if
        
    end subroutine

    subroutine enable_warnings(value)
        use m_common, only: disable_warnings
        implicit none
        logical, intent(in) :: value
        disable_warnings = .not. value
    end subroutine
#endif

end module m_sw_mono
