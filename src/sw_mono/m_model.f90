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
!> @file m_sw_mono.f90
!> @brief This file includes the m_model module.

!> @brief Module m_model.
module m_model
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_obs
!     use m_time_screen                                                                                              !NOADJ
    implicit none
#ifndef CPP_ADJ
    include "mpif.h"
    include 'dmumps_struc.h'
#endif
   
   !===================================================================================================================!
   !  Discrete Model Unknows Structure
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
    
    type Model

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
        !> Froude threshold for LPI method
        real(rp) :: frLPI
        !> Gravity acceleration (m2/s)
        real(rp) :: gravity
        !> Depth threshold
        real(rp) :: heps
        !> Parameter for LPI method
        real(rp) :: mLPI
        !> Implicit coefficient for Preissmann scheme
        real(rp) :: theta_preissmann
        !> Bounday conditions
        type(BoundaryCondition), dimension(:), allocatable :: bc
	!>scale for the implicit diffusive wave
	character(len=16) :: scale_idw
      
    end type


#ifndef CPP_ADJ
    interface
!   
!         subroutine apply_control(mdl, msh, dof, ctrl)
!             import Control
!             import Mesh
!             import Model
!             import Unknowns
!             implicit none
!             type(Control), intent(inout) :: ctrl
!             type(Model), intent(inout) :: mdl
!             type(Mesh), intent(inout) :: msh
!             type(Unknowns), intent(inout) :: dof
!         end subroutine apply_control
!         subroutine calc_cost(mdl, msh, imp, dof, ctrl, obs, cost)
!             use m_control
!             import Control
!             import ImplicitMatrix
!             import Mesh
!             import Model
!             import Observations
!             import rp
!             import Unknowns
!             implicit none
!             type(Model), intent(inout) :: mdl
!             type(Mesh), intent(inout) :: msh
!             type(ImplicitMatrix), intent(inout) :: imp
!             type(Unknowns), intent(inout) :: dof
!             type(Control), intent(inout) :: ctrl
!             type(Observations), optional :: obs
!             real(rp), intent(out) :: cost
!         end subroutine
        subroutine generate_observations(mdl, msh, imp, dof, obs)
            import ImplicitMatrix
            import Mesh
            import Model
            import Observations
            import Unknowns
            implicit none
            type(Model), intent(inout) :: mdl
            type(Mesh), intent(inout) :: msh
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
            type(Observations), intent(inout), optional :: obs
        end subroutine
        subroutine preissmann_timestep(mdl, msh, imp, dof)
            import ImplicitMatrix
            import Mesh
            import Model
            import Unknowns
            implicit none
            type(Model), intent(in) :: mdl
            type(Mesh), intent(inout) :: msh
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
        end subroutine
	subroutine implicit_diffusive_wave(mdl, msh, imp, dof)
            import ImplicitMatrix
            import Mesh
            import Model
            import Unknowns
            implicit none
            type(Model), intent(in) :: mdl
            type(Mesh), intent(inout) :: msh
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
        end subroutine
        subroutine time_loop(mdl, msh, imp, dof, obs)
            import ImplicitMatrix
            import Mesh
            import Model
            import Observations
            import Unknowns
            implicit none
            type(Model), intent(inout) :: mdl
            type(Mesh), intent(inout) :: msh
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
            type(Observations), intent(inout), optional :: obs
        end subroutine
    end interface
#endif

    
! ! ! 
! ! ! 
! ! ! 
! ! !    real(rp) bathy_slopemean                                    !< mean slope of the canal riverbed
! ! !    real(rp), dimension(:), allocatable  ::  bathy_cell         !< Bathymetry at msh cells gravity center.
! ! !    integer(ip) :: nb_bathy_control_pts
! ! ! 
! ! !    real(rp), dimension(:), allocatable  ::  tflux1             !< Flux1 for finite volume scheme.
! ! ! 
! ! !    real(rp), dimension(:), allocatable  ::  tflux2             !< Flux2 for finite volume scheme.
! ! ! 
! ! !    real(rp), dimension(:), allocatable  ::  surfaceWidth       !< Surface width.
! ! ! 
! ! !    real(rp), dimension(:), allocatable  :: hn, hc                  !< Normal depth, critical depth
! ! !    real(rp), dimension(:), allocatable  :: backwater                 !TMP
! ! !    
! ! !    !===================================================================================================================!
! ! !    !  Boundary Condition Structures
! ! !    !===================================================================================================================!
! ! ! 
! ! !    !> Structure of timeseries.
! ! !    !!
! ! !    !! \details This structure is used for the definition of timeseries (i.e for each \f$t\f$, a \f$Y(t)\f$ is
! ! !    !! given.
! ! !    type timeseries
! ! !       real(rp), dimension(:), allocatable  ::  t  !< Time array.
! ! !       real(rp), dimension(:), allocatable  ::  y  !< Flow array.
! ! !    end type
! ! ! 
! ! !    !> Structure of hydrograph.
! ! !    !!
! ! !    !! \details This structure is used for the definition of temporal hydrograph (i.e for each \f$t\f$, a \f$Q_{in}(t)\f$ is
! ! !    !! given.
! ! !    type hydrograph
! ! !       real(rp), dimension(:), allocatable  ::  t  !< Time array.
! ! !       real(rp), dimension(:), allocatable  ::  q  !< Flow array.
! ! !    end type
! ! ! 
! ! !    !> Structure of Fourier series hydrograph.
! ! !    !!
! ! !    !! \details This structure is used for the definition of Fourier series hydrograph. Hydrograph is defined thus that
! ! !    !!
! ! !    !! \f$Q_{in}(t)=\frac{a_{0}}{2}+\sum_{n=1}^{N_{FS}}\left(a_{n}\cos(nt\frac{2\pi}{T})+b_{n}\sin(nt\frac{2\pi}{T})\right)\f$ 
! ! !    !! with \f$\{a_{0};a_{n},b_{n}\} \f$ , \f$n\in[1..N_{FS}]\f$, are the Fourier coefficients and \f$T\f$ is the total 
! ! !    !! simulation time. The lower frequency represented by the Fourier serie is \f$1/T\f$ and the highest one 
! ! !    !! is \f$N_{FS}/T\f$. 
! ! !    type hydrograph_FS
! ! !    
! ! !       real(rp) :: a0                                  !< Variable \f$a_{0}\f$.
! ! !       real(rp), dimension(:), allocatable  ::  a      !< Variables \f$a_{i}\f$.
! ! !       real(rp), dimension(:), allocatable  ::  b      !< Variables \f$b_{i}\f$.
! ! !       integer(ip), dimension(:), allocatable  ::  n   !< Variable \f$N_{FS}\f$.
! ! !    end type
! ! ! 
! ! !    !> Structure of ratcurve
! ! !    !!
! ! !    !! \details This structure is used for the definition of ratcurve: \f$h(q)\f$.
! ! !    type ratcurve
! ! !       real(rp), dimension(:), allocatable  ::  h !< Height variable.
! ! !       real(rp), dimension(:), allocatable  ::  q !< Flow variable.
! ! !    end type
! ! ! 
! ! !    type hydrograph_qlat
! ! !       real(rp), dimension(:), allocatable  ::  t !< Time as read from text files
! ! !       real(rp), dimension(:), allocatable  ::  tsim !< Discretized time
! ! !       real(rp), dimension(:,:), allocatable  ::  q !< Lateral inflows as read from text files
! ! !       real(rp), dimension(:), allocatable  ::  qsim !< Lateral inflows discretized in time
! ! !       real(rp), dimension(:), allocatable  ::  x !< Abscisse of lateral inflows
! ! !       integer(ip), dimension(:), allocatable  ::  loc !< Mesh index of lateral inflows
! ! ! !       character(len=25), dimension(:), allocatable  ::  id !< Name of text files
! ! !       real(rp), dimension(:,:), allocatable  ::  params !< Inversion parameters : sigma, mu, alpha_regul, read from idlist
! ! !    end type   
! ! ! 
! ! !    !> Structure of boundaries conditions
! ! !    !!
! ! !    !! \details This structure is used for the definition of boundaries conditions (type of hydrograph, of ratcurve, ...).
! ! !    type bcs
! ! ! 
! ! !       character(len=lchar) ::  typehyd !< Type of input hydrograph used: 'file_FS' for Fourier series definition
! ! !    !! or 'file' for temporal definition.
! ! !       
! ! !       character(len=lchar) ::  typerat !< Type of ratcurve used: 'file' for a definition \f$H(Q)\f$, else definition with
! ! !    !! \f$\alpha h ^{\beta}\f$.
! ! ! 
! ! !       type(hydrograph) ::  hyd         !< Temporal hydrograph type.
! ! ! 
! ! !       type(timeseries) ::  dsz         !< Time series of downstream elevations
! ! ! 
! ! !       type(hydrograph_FS) ::  hyd_FS   !< Fourier Series hydrograph type.
! ! ! 
! ! !       type(ratcurve) ::  rat           !< Ratcurve type.
! ! ! 
! ! !       type(hydrograph_qlat) ::  hyd_lat         !< Temporal hydrograph type for lateral flow.
! ! !       
! ! !    end type bcs
! ! ! 
! ! !    type( bcs )  ::  bc                 !< Boundaries conditions definitions.
! ! ! 
! ! !    !===================================================================================================================!
! ! !    !  Strickler (Manning in the code) coeff model parameters
! ! !    !===================================================================================================================!
! ! !    !> Structure of Strickler coeff
! ! !    !!
! ! !    !! \details This structure is used for the definition of K(x), the Strickler coefficient modeling the friction
! ! !    type strickler_params
! ! !       integer(ip) :: nb_diff_K_reachs       !< nb of reachs along which the strickler params are constant
! ! !       integer(ip) :: nb_K_params            !< nb of parameters defining the Strickler coeff K
! ! !       integer(ip), dimension(:), allocatable :: cs_start, cs_end   !< start and end cross-sections of each separate reach with different Strickler parameters
! ! !       real(rp), dimension(:), allocatable :: alpha, beta        !< TEST : for now, the model is K(h,x)=alpha(x)*h^beta(x)
! ! !       real(rp), dimension(:, :), allocatable :: coeffs        !< TEST : for now, the model is K(h,x)=alpha(x)*h^beta(x)
! ! !    end type strickler_params
! ! ! 
! ! !    type( strickler_params ) :: K_params
! ! ! 
! ! !    !==============================================================================================!
! ! !    !  Bathymetry points (for piecewise models)
! ! !    !==============================================================================================!
! ! !    integer(ip), dimension(:), allocatable  :: bathy_indices
! ! !    real(rp), dimension(:), allocatable  :: bathy_points
! ! ! 
! ! !    !==============================================================================================!
! ! !    !  First guess of model parameters used as control in the DA process
! ! !    !==============================================================================================!
! ! !    real(rp), dimension(:), allocatable  :: bathy_cell_fg
! ! !    real(rp), dimension(:), allocatable  :: bathy_points_fg
! ! !    real(rp), dimension(:), allocatable :: alpha_K_fg
! ! !    real(rp), dimension(:), allocatable :: beta_K_fg
! ! !    real(rp), dimension(:), allocatable  :: qin_fg
! ! !    real(rp), dimension(:,:), allocatable  :: qlat_fg
! ! !    
! ! !    !==============================================================================================!
! ! !    !  Variables v, when the variable change v = Bdemi*(u-u_b) is used
! ! !    !==============================================================================================!
! ! !    ! /!\ For now, there are no cross-covariances (Bdemi is only block diagonal). 
! ! !    real(rp), dimension(:), allocatable  :: bathy_cell_chg
! ! !    real(rp), dimension(:), allocatable  :: bathy_points_chg
! ! !    real(rp), dimension(:), allocatable :: alpha_K_chg
! ! !    real(rp), dimension(:), allocatable :: beta_K_chg
! ! !    real(rp), dimension(:), allocatable  :: qin_chg 
! ! !    real(rp), dimension(:,:), allocatable  :: qlat_chg
! ! ! 
! ! !    !==============================================================================================!
! ! !    !  Covariances matrices
! ! !    !==============================================================================================!
! ! ! 
! ! !    ! Bdemi -> /!\ contain lower and upper triangular matrices from Cholesky decomposition
! ! !    ! So Bdemi*Bdemi != B, need to select only lower or upper half
! ! !       real(rp), dimension(:, :), allocatable  ::  Bdemi_hydrograph !< Block of \f$ B^{1/2} \f$ for the hydrograph.
! ! !       real(rp), dimension(:, :, :), allocatable  ::  Bdemi_hydrograph_lat !< Block of \f$ B^{1/2} \f$ for the lateral hydrograph.
! ! !    real(rp), dimension(:, :), allocatable  ::  Bdemi_bathy      !< Block of \f$ B^{1/2} \f$ for the bathymetry.
! ! ! !!! TO DO, PAG clean this doublon
! ! !    real(rp), dimension(:), allocatable  ::  Bdemi_K_alpha    !< Block of \f$ B^{1/2} \f$ for alpha (just an array because we assume Bdemi is diagonal)
! ! !    real(rp), dimension(:), allocatable  ::  Bdemi_K_beta     !< Block of \f$ B^{1/2} \f$ for beta (just an array because we assume Bdemi is diagonal)
! ! !    real(rp), dimension(:,:), allocatable  ::  Bdemi_K_alpha_m    !< Block of \f$ B^{1/2} \f$ for alpha 
! ! !    real(rp), dimension(:,:), allocatable  ::  Bdemi_K_beta_m  !< Block of \f$ B^{1/2} \f$ for beta 
! ! ! 
! ! ! 
! ! ! 
! ! !    ! Bdemi inv -> same remark as for Bdemi
! ! !    real(rp), dimension(:, :), allocatable  ::  Bdemi_hydro_inv !< Block of \f$ B^{-1/2} \f$ for the hydrograph.
! ! !    real(rp), dimension(:, :), allocatable  ::  Bdemi_bathy_inv      !< Block of \f$ B^{-1/2} \f$ for the bathymetry.
! ! !    
! ! !    !===================================================================================================================!
! ! !    !  Recording Structures
! ! !    !===================================================================================================================!
! ! !    !> Structure of stations.
! ! !    !!
! ! !    !! \details This structure is used for the definition of observation station. In this structure are defined lenght
! ! !    !! of the observation station, time step, weight, index observation time, ...
! ! ! 
! ! !    type station_obs
! ! ! 
! ! !       type( point_in_mesh ), dimension(:), allocatable  ::  pt    !< Array of point in the msh.
! ! ! 
! ! !       real(rp)  :: weight                                !< Weight of observations.
! ! ! 
! ! !       real(rp)  :: length                                !< Length of river.
! ! !       
! ! !       real(rp)  :: dt_offset                             !< Time of first observation (or offset).
! ! ! 
! ! !       real(rp)  :: dt                                    !< Frequency of observation ( satellite time repetitiveness).
! ! ! 
! ! !       real(rp),dimension(:), allocatable :: dt_obs       !< Array with observation time.
! ! ! 
! ! !       integer(ip) :: ind_t                               !< Index observation time.
! ! ! 
! ! !       integer(ip) :: nb_dt                               !< Number of observation time.
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  t         !< Time observation \f$T^{obs}\f$. 
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  h         !< Water height observed: \f$\bar{h}^{obs}\f$.
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  s         !< Surface observed: \f$\bar{S}^{obs}\f$.
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  q         !< Flow  observed: \f$\bar{Q}^{obs}\f$.
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  h_bathy   !< Water elevation observed: \f$\bar{H}^{obs}\f$.
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  w         !< surface width observed: \f$\bar{w}^{obs}\f$.
! ! ! 
! ! !       real(rp), dimension(:), allocatable  ::  slope     !< surface slope observed
! ! ! 
! ! ! 
! ! !    end type station_obs
! ! ! 
! ! !    !===================================================================================================================!
! ! !    !  Recording Variables
! ! !    !===================================================================================================================!
! ! ! 
! ! !    type( station_obs ), dimension(:), allocatable  ::  station !< Observation station definition.
! ! ! 
! ! ! 
! ! ! 
! ! !    !===================================================================================================================!
! ! !    !  Input variables specific to model (in addition to m_common)
! ! !    !===================================================================================================================!
! ! !    
! ! !    character(len=lchar) :: crosssection_extrapolation !< Model for extrapolating above highest triplet (Z,W,Y)
! ! ! 
! ! !    real(rp)     ::  g                                 !< Gravity constant \f$g\f$.
! ! !    real(rp)     ::  heps                              !< Minimal water depth.
! ! !    real(rp)     ::  qeps                              !< Minimal discharge.
! ! !    real(rp)     ::  theta_preissman                   !< Implicit factor for Preissmann scheme
! ! !    integer(ip)  ::  friction                          !< Activation of a Friction Law in Model.
! ! !    character(len=lchar)  :: K_model                   !< Choice of law for the Strickler coeff K
! ! ! 
! ! !    integer(ip)  ::  c_manning                         !< Activation of Manning \f$K\f$ (if ==1) in the control vector 
! ! ! !!\f$\mathbf{k}\f$. 
! ! !    integer(ip)  ::  c_bathy                           !< Activation of bathymetry \f$z_b\f$ (if ==1) in the control 
! ! ! !!vector \f$\mathbf{k}\f$. 
! ! !    integer(ip)  ::  c_ic                              !< Activation of initiales conditions \f$S^0\f$ and \f$Q^0\f$ 
! ! ! !! (if ==1) in the control vector \f$\mathbf{k}\f$. 
! ! !    integer(ip)  ::  c_hydrograph                      !< Activation of hydrograph \f$Q_{in}(t)\f$ in the control vector 
! ! ! !! \f$\mathbf{k}\f$. 
! ! ! !! - if ==1: Using of temporal form \f$Q_{in}(t^{i})\f$.
! ! ! !! - if ==2: Using of Fourier series form \f$Q_{in}(t)=\frac{a_{0}}{2}+\sum_{n=1}^{N_{FS}}.
! ! ! !! \left(a_{n}\cos(nt\frac{2\pi}{T})+b_{n}\sin(nt\frac{2\pi}{T})\right)\f$.
! ! ! 
! ! !    integer(ip)  ::  c_ratcurve                        !< Activation of ratcurve \f$H_{out}(Q)\f$ in the control vector 
! ! ! !!\f$\mathbf{k}\f$. 
! ! ! !! - if ==1: Using of the tradional ratcurve \f$H_{out}=f(Q)\f$.
! ! ! !! - if ==2: Using of the relation \f$H_{out}=\alpha Q ^{\beta}\f$.
! ! ! 
! ! !    real(rp)     ::  alpha_lbound                      !< Lower bound for the Manning \f$\alpha\f$ \f$K=\alpha \dot h^{\beta}\f$.
! ! !    real(rp)     ::  alpha_ubound                      !< Upper bound for the Manning \f$\alpha\f$ \f$K=\alpha \dot h^{\beta}\f$.
! ! !    real(rp)     ::  beta_lbound                       !< Lower bound for the Manning \f$\beta\f$ \f$K=\alpha \dot h^{\beta}\f$.
! ! !    real(rp)     ::  beta_ubound                       !< Upper bound for the Manning \f$\beta\f$ \f$K=\alpha \dot h^{\beta}\f$.
! ! !    real(rp)     ::  bathy_lbound                      !< Lower bound for the bathymetry \f$z_b\f$.
! ! !    real(rp)     ::  bathy_ubound                      !< Upper bound for the bathymetry \f$z_b\f$.
! ! !    real(rp)     ::  bathy_delta                       !< Bounds interval around initial values for the bathymetry \f$z_b\f$.
! ! !    real(rp)     ::  hydrograph_lbound                 !< Lower bound for the hydrograph \f$Q_{in}(t)\f$.
! ! !    real(rp)     ::  hydrograph_ubound                 !< Upper bound for the hydrograph \f$Q_{in}(t)\f$.
! ! !    real(rp)     ::  alpha_ratcurve_lbound             !< Lower bound for the rating curve \f$\alpha\f$
! ! !    real(rp)     ::  alpha_ratcurve_ubound             !< Upper bound for the rating curve \f$\alpha\f$
! ! !    real(rp)     ::  beta_ratcurve_lbound              !< Lower bound for the rating curve \f$\alpha\f$
! ! !    real(rp)     ::  beta_ratcurve_ubound              !< Upper bound for the rating curve \f$\alpha\f$
! ! ! 
! ! !    real(rp)     ::  eps_manning                       !< Value of \f$\epsilon\f$ in perturbation control vector 
! ! ! !!\f$J(\mathbf{k}+\epsilon \delta \mathbf{k})\f$ for the Manning \f$K\f$.
! ! !    real(rp)     ::  eps_bathy                         !< Value of \f$\epsilon\f$ in perturbation control vector 
! ! ! !!\f$J(\mathbf{k}+\epsilon \delta \mathbf{k})\f$ for the bathymetry \f$z_b\f$.
! ! !    real(rp)     ::  eps_ic                            !< Value of \f$\epsilon\f$ in perturbation control vector 
! ! ! !!\f$J(\mathbf{k}+\epsilon \delta \mathbf{k})\f$ for the initiales conditions \f$S^0\f$ and \f$Q^0\f$.
! ! !    real(rp)     ::  eps_hydrograph                    !< Value of \f$\epsilon\f$ in perturbation control vector 
! ! ! !!\f$J(\mathbf{k}+\epsilon \delta \mathbf{k})\f$ for the  hydrograph \f$Q_{in}(t)\f$.
! ! !    real(rp)     ::  eps_ratcurve                      !< Value of \f$\epsilon\f$ in perturbation control vector 
! ! ! !!\f$J(\mathbf{k}+\epsilon \delta \mathbf{k})\f$ for the ratcurve \f$H_{out}(Q)\f$.
! ! ! 
! ! !    character(len=lchar)  :: bathy_regul_type          !< Type of regularization for the bathymetry
! ! !    character(len=lchar)  :: strickler_regul_type      !< Type of regularization for the Strickler coefficient
! ! !    character(len=lchar)  :: hydro_regul_type          !< Type of regularization for the hydrogram
! ! !    character(len=lchar)  :: K_covariance              !< Type of regularization for K
! ! ! 
! ! ! 
! ! !    real(rp)     ::  regul_manning                     !< Weight of the Manning regularization term in the cost function.
! ! !    real(rp)     ::  regul_bathy                       !< Weight of the bathymetry regularization term in the 
! ! ! ! cost function.
! ! !    real(rp)     ::  regul_ic                          !< Weight of the initiales conditions regularization term in the 
! ! ! ! cost function.
! ! !    real(rp)     ::  regul_hydrograph, regul_hydrographlat                  !< Weight of the hydrograph regularization term in the 
! ! ! ! cost function.
! ! !    real(rp)     ::  regul_ratcurve                    !< Weight of the ratcurve regularization term in the
! ! !    real(rp)     ::  regul_gamma                      !< Factor for the computation of the iterative regularization alpha
! ! !    
! ! !    logical  ::  var_chg               !< if true, does the change of variables u = u_fg + Bdemi*v,
! ! !    logical :: norm_grad0
! ! !    ! where u_fg is the first guess, Bdemi is the square root of the background error cov matrix,
! ! !    ! u the model parameters and v the new variable
! ! !    
! ! !    integer(ip) :: run_status                          !< 0 is run is sucessful, negative otherwise
! ! ! 
! ! !    integer(ip)  ::  fix_time_step_serie               !< Variable to fix the time step serie.
! ! ! 
! ! !    character(len=lchar)  ::  bathymetry               !< Type of bathymetry model ('parametric', 'piecewise_constant', 'piecewise_linear')
! ! !    integer(ip)  ::  bathy_first, bathy_last, bathy_step
! ! !    real(rp)     ::  sigma_bathy, mu_bathy             !< Standard deviation and correlation length for the bathymetry.
! ! !    real(rp)     ::  sigma_alpha                       !< Standard deviation for the alpha manning.
! ! !    real(rp)     ::  sigma_beta                        !< Standard deviation for the beta manning.
! ! !    real(rp)     ::  mu_alpha,mu_beta                  !< correlation length for for the alpha manning and beta manning 
! ! !    real(rp)     ::  sigma_hydrograph, mu_hydrograph   !< Standard deviation and correlation length for the hydrograph.
! ! !    real(rp)     ::  sigma_hydrographlat, mu_hydrographlat   !< Standard deviation and correlation length for the hydrograph.
! ! !    
! ! !    
! ! !    real(rp)     ::  backwater_hmin, backwater_alpha
! ! !    real(rp)     ::  normal_depth_heps
! ! ! 
! ! !    integer(ip)  ::  ite_period_alpha_reg_calc !< the minimization algo recomputes the regularization coeff alpha every ite_period_alpha_reg_calc iterations - if = 0 does not recompute alpha
! ! !    
! ! !    integer(ip)  ::  zb_cs_1, zb_cs_2
! ! !    
! ! !    ! BACKWARDS COMPATIBILITY
! ! !    real(rp) :: alpha, beta                             !< Constant in space Strikler terms
! ! !    character(len=lchar)  ::  regul_type               !< Type of regulation ('tikhonov_h0' only so far)
! ! !    
! ! !    
! ! !    integer :: debug_flag
! ! !    
! ! !    !===================================================================================================================!
! ! !    !  Input variables namelist (m_common + model specific upper ones)
! ! !    !===================================================================================================================!
! ! ! 
! ! !    namelist/list_input/ &
! ! !       mesh_name,&
! ! !       crosssection_extrapolation,&
! ! !       remesh,&
! ! !       smooth,&
! ! !       addPointsW,&
! ! !       addPointsE,&
! ! !       post,&
! ! !       bc_W,&
! ! !       bc_E,&
! ! !       lat_inflow,&
! ! !       initialization, &
! ! !       backwater_hmin, & ! DEPRECATED ?
! ! !       backwater_alpha, & ! DEPRECATED ?
! ! !       normal_depth_heps, &
! ! !       bathymetry, &
! ! !       bathy_first, &
! ! !       bathy_last, &
! ! !       bathy_step, &
! ! !       w_exact,&
! ! !       ts,&
! ! !       dtw,&
! ! !       dtp,&
! ! !       dta,&
! ! !       spatial_scheme, &
! ! !       theta_preissman, &
! ! !       dt,&
! ! !       adapt_dt,&
! ! !       cfl,&
! ! !       friction,&
! ! !       K_model,&
! ! !       alpha_ratcurve,&
! ! !       beta_ratcurve,&
! ! !       w_froude,&
! ! !       steady_state,&
! ! !       eps_q_steady_state,&
! ! !       eps_s_steady_state,&
! ! !       g,&
! ! !       heps, &
! ! !       qeps, &
! ! !       Qinterface, &
! ! !       w_obs,&
! ! !       use_obs,&
! ! !       eps_min,&
! ! !       eps_morozov,&
! ! !       c_manning, &
! ! !       c_bathy, &
! ! !       c_ic, &
! ! !       c_hydrograph, &
! ! !       c_ratcurve, &
! ! !       
! ! ! 
! ! !       zb_cs_1, &
! ! !       zb_cs_2, &
! ! !       
! ! !       alpha_lbound,&
! ! !       alpha_ubound,&
! ! !       beta_lbound,&
! ! !       beta_ubound,&
! ! !       bathy_lbound,&
! ! !       bathy_ubound,&
! ! !       bathy_delta,&
! ! !       hydrograph_lbound,&
! ! !       hydrograph_ubound,&
! ! !       alpha_ratcurve_lbound,&
! ! !       alpha_ratcurve_ubound,&
! ! !       beta_ratcurve_lbound,&
! ! !       beta_ratcurve_ubound,&
! ! ! 
! ! !       eps_manning, &
! ! !       eps_bathy, &
! ! !       eps_ic, &
! ! !       eps_hydrograph, &
! ! !       eps_ratcurve, &
! ! ! 
! ! !       regul_manning, &
! ! !       regul_bathy, &
! ! !       regul_ic, &
! ! !       regul_hydrograph, &
! ! !       regul_ratcurve, &
! ! !       var_chg, &
! ! !       norm_grad0, &
! ! ! 
! ! !       strickler_regul_type,&
! ! !       bathy_regul_type, &
! ! !       hydro_regul_type, &
! ! !       K_covariance,&
! ! !       
! ! !       sigma_bathy, &
! ! !       mu_bathy, &
! ! !       sigma_alpha, &
! ! !       sigma_beta, &
! ! !       mu_alpha, &
! ! !       mu_beta, &
! ! !       sigma_hydrograph, &
! ! !       mu_hydrograph, &
! ! !       sigma_hydrographlat, &
! ! !       mu_hydrographlat, &
! ! ! 
! ! !       ite_period_alpha_reg_calc, &
! ! !       regul_gamma, &
! ! !       
! ! !       ! BACKWARDS COMPATIBILITY
! ! !       alpha, &
! ! !       beta, &
! ! !       regul_type
      
#ifndef CPP_ADJ
    contains
    
    
    subroutine set_timeseries(bc, t, y)
        implicit none
        type(BoundaryCondition), intent(inout) :: bc
        real(rp), dimension(:), intent(in) :: t
        real(rp), dimension(:), intent(in) :: y

        if (size(t) /= size(y)) then
            call f90wrap_abort("t and y must be of same size")
        end if
        
        if (allocated(bc%ts%t)) deallocate(bc%ts%t)
        if (allocated(bc%ts%y)) deallocate(bc%ts%y)
        
        allocate(bc%ts%t(size(t)))
        allocate(bc%ts%y(size(t)))
        bc%ts%t(:) = t
        bc%ts%y(:) = y
        
    end subroutine
    
    
    subroutine implicitmatrix_initialise(imp, mdl, msh)
        implicit none
        type(Model), intent(in) :: mdl
        type(Mesh), intent(in) :: msh
        type(ImplicitMatrix), intent(out) :: imp
        
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
    
    
    subroutine model_initialise(mdl, msh)
        implicit none
        type(Mesh), intent(in) :: msh
        type(Model), intent(out) :: mdl
        
        ! Index of boundary condition
        integer(ip) :: ibc
        ! Index of segment
        integer(ip) :: iseg
        ! Number of boundary conditions
        integer(ip) :: nbc
        
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
        allocate(mdl%bc(nbc))
        do ibc = 1, nbc
            mdl%bc(ibc)%id = "undefined"
        end do
        
        ! Set default values
        mdl%tc = 0.0
        mdl%ts = 0.0
        mdl%te = 1.0
        mdl%dt = 1.0
        mdl%dtout = 1.0
        mdl%heps = 0.01
        mdl%frLPI = 0.7
        mdl%mLPI = 10.0
        mdl%gravity = 9.81
        mdl%theta_preissmann = 0.7
	mdl%scale_idw = "large"
        
    end subroutine
    
    
    subroutine model_finalise(mdl)
        implicit none
        type(Model), intent(inout) :: mdl
        
        ! Deallocate boundary conditions
        if (allocated(mdl%bc)) then
            ! TODO finalise bcs
            deallocate(mdl%bc)
        end if
        
    end subroutine
    
    
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
#endif

end module m_model
