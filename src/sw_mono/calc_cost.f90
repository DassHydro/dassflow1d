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
subroutine calc_cost(cost, mdl, ctrl, obs)
    use m_common
    use m_control, only: Control
    use m_mesh, only: Mesh
    use m_sw_mono, only: ImplicitMatrix
    use m_sw_mono, only: Model
    use m_sw_mono, only: Unknowns
    use m_numeric
    use m_obs, only: Observations
!     use m_obs, only: cost_L2_heights
    use m_obs, only: cost_R2_heights
#ifndef CPP_ADJ
    use m_obs, only: data_from_stations
    use m_obs, only: write_observations_status
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param mdl Model
    type(Model), intent(inout) :: mdl
    !> @param ctrl Control
    type(Control), intent(inout) :: ctrl
    !> @param obs Observations
    type(Observations), intent(inout) :: obs
    !> @param cost Cost function
    real(rp), intent(out) :: cost
    !=LOCAL VARIABLES==================================================================================================!
    integer(ip) :: ics!,icsm,icsp
    !integer(ip) :: iseg
#ifdef REGUL_D2_ALPHA
#ifndef REGUL_ALPHA
#define REGUL_ALPHA
#define REGUL_GRAD
#ifndef REGUL
#define REGUL
#endif
#endif
#endif    
#ifdef REGUL_D1_ALPHA
#ifndef REGUL_ALPHA
#define REGUL_ALPHA
#define REGUL_GRAD
#ifndef REGUL
#define REGUL
#endif
#endif
#endif    
#ifdef REGUL_ALPHA
!     integer(ip) :: ics
    real(rp) :: alphaL, alphaC, alphaR
#endif    
#ifdef REGUL_BATHY
#ifndef REGUL
#define REGUL
#endif
#define REGUL_GRAD
!     integer(ip) :: ics
    real(rp) :: bL, bC, bR
#endif    
    real(rp) :: sumreg
!     print *, "CALC_COST!"
#ifdef REGUL_GRAD
    real(rp) :: gradL, gradR
#endif    
    
#ifndef CPP_ADJ
    ! CHECK-UP
    if (mdl%scheme == "undefined") call abort_solver("Scheme is undefined")
#endif
    
#ifndef CPP_ADJ
    ! Initialise obs/est arrays if necessary
    if (.not.allocated(obs%obs)) then
        if (.not.allocated(obs%stations)) then
            call abort_solver("No observations defined")
        end if
        call data_from_stations(obs, mdl%discharge_estimation)
    end if
    
#endif


    ! Apply control
    call apply_control(ctrl, mdl)
    
    ! Apply fields
!     if (allocated(mdl%msh%bathy_field%x)) then
    call apply_bathy_field(mdl%msh, 1.1_rp * mdl%heps)
!     end if
!     if (allocated(mdl%msh%seg(mdl%msh%nseg)%strickler_fields)) then
    call apply_strickler_fields(mdl%msh)
!     end if

    if (mdl%steady_states) then
    
        ! TODO ensure all timeseries (BCs, inflow, obs) have same time occurences
        
        ! Steady states loop
        call steady_states_loop(mdl, obs)
    
    else
    
        ! Time loop
        call time_loop(mdl, obs)
        
    end if

    ! Compute cost
#ifndef CPP_ADJ
    if (.not. allocated(obs%w)) then
        call abort_solver("Observations weights not set")
    end if
#endif
!     cost = cost_L2_heights(obs%est(1, :), obs%obs(1, :))
    cost = cost_R2_heights(obs%est(1, :), obs%obs(1, :), obs%w(1, :))
    
#ifndef CPP_ADJ
    if (mdl%observations_status) then
        call write_observations_status(obs)
    end if
#endif
    
#ifdef REGUL_D2_ALPHA
#ifndef CPP_ADJ
    if (mdl%msh%nseg > 1) call abort_solver("Regul(d2dx(alpha)) only works on single segment domain")
#endif
    ics = mdl%msh%seg(1)%first_cs
    alphaR = mdl%msh%cs(ics+1)%strickler_params(1)
    alphaC = mdl%msh%cs(ics)%strickler_params(1)
    gradR = (alphaR - alphaC) / mdl%msh%cs(ics+1)%deltademi
    gradL = 0.0
    sumreg = ((gradR - gradL) / mdl%msh%cs(ics)%delta)**2
    do ics = mdl%msh%seg(1)%first_cs+1, mdl%msh%seg(1)%last_cs-1
        alphaR = mdl%msh%cs(ics+1)%strickler_params(1)
        alphaC = mdl%msh%cs(ics)%strickler_params(1)
        alphaL = mdl%msh%cs(ics-1)%strickler_params(1)
        gradR = (alphaR - alphaC) / mdl%msh%cs(ics+1)%deltademi
        gradL = (alphaC - alphaL) / mdl%msh%cs(ics)%deltademi
        sumreg = sumreg + ((gradR - gradL) / mdl%msh%cs(ics)%delta)**2
    end do
    ics = mdl%msh%seg(1)%last_cs
    alphaC = mdl%msh%cs(ics)%strickler_params(1)
    alphaL = mdl%msh%cs(ics-1)%strickler_params(1)
    gradR = 0.0
    gradL = (alphaC - alphaL) / mdl%msh%cs(ics)%deltademi
    sumreg = sumreg + ((gradR-gradL) / mdl%msh%cs(ics)%delta)**2

#endif

#ifdef REGUL_D1_ALPHA
#ifndef CPP_ADJ
    if (mdl%msh%nseg > 1) call abort_solver("Regul(d1dx(alpha)) only works on single segment domain")
#endif
    ics = mdl%msh%seg(1)%first_cs
    alphaR = mdl%msh%cs(ics+1)%strickler_params(1)
    alphaC = mdl%msh%cs(ics)%strickler_params(1)
    gradR = (alphaR - alphaC) / mdl%msh%cs(ics+1)%deltademi
    gradL = 0.0
    sumreg = (0.5 * (gradR + gradL))**2
    do ics = mdl%msh%seg(1)%first_cs+1, mdl%msh%seg(1)%last_cs-1
        alphaR = mdl%msh%cs(ics+1)%strickler_params(1)
        alphaC = mdl%msh%cs(ics)%strickler_params(1)
        alphaL = mdl%msh%cs(ics-1)%strickler_params(1)
        gradR = (alphaR - alphaC) / mdl%msh%cs(ics+1)%deltademi
        gradL = (alphaC - alphaL) / mdl%msh%cs(ics)%deltademi
        sumreg = sumreg + (0.5 * (gradR + gradL))**2
    end do
    ics = mdl%msh%seg(1)%last_cs
    alphaC = mdl%msh%cs(ics)%strickler_params(1)
    alphaL = mdl%msh%cs(ics-1)%strickler_params(1)
    gradR = 0.0
    gradL = (alphaC - alphaL) / mdl%msh%cs(ics)%deltademi
    sumreg = sumreg + (0.5 * (gradR + gradL))**2
#endif


#ifdef REGUL_D1_BATHY
#ifndef CPP_ADJ
    if (mdl%msh%nseg > 1) call abort_solver("Regul(d1dx(bathy)) only works on single segment domain")
#endif
    ics = mdl%msh%seg(1)%first_cs
    bR = mdl%msh%cs(ics+1)%bathy
    bC = mdl%msh%cs(ics)%bathy
    gradR = (bR - bC) / mdl%msh%cs(ics+1)%deltademi
    gradL = 0.0
    sumreg = (0.5 * (gradR + gradL))**2
    do ics = mdl%msh%seg(1)%first_cs+1, mdl%msh%seg(1)%last_cs-1
        bR = mdl%msh%cs(ics+1)%bathy
        bC = mdl%msh%cs(ics)%bathy
        bL = mdl%msh%cs(ics-1)%bathy
        gradR = (bR - bC) / mdl%msh%cs(ics+1)%deltademi
        gradL = (bC - bL) / mdl%msh%cs(ics)%deltademi
        sumreg = sumreg + (0.5 * (gradR + gradL))**2
    end do
    ics = mdl%msh%seg(1)%last_cs
    bC = mdl%msh%cs(ics)%bathy
    bL = mdl%msh%cs(ics-1)%bathy
    gradR = 0.0
    gradL = (bC - bL) / mdl%msh%cs(ics)%deltademi
    sumreg = sumreg + (0.5 * (gradR + gradL))**2
#endif

    mdl%cost_obs = cost
#ifdef REGUL
    if (mdl%gamma_reg < 0.0_rp) then
        if (sumreg > 1e-12) then
            mdl%gamma_reg = 0.1 * cost / sumreg
#ifndef CPP_ADJ
            if (.not. mdl%disable_stdout) then 
                print *, "INITIALISATION gamma_reg=", mdl%gamma_reg
            end if
#endif
            mdl%cost_reg = mdl%gamma_reg * sumreg
        else
#ifndef CPP_ADJ
            if (.not. mdl%disable_stdout) then 
                print *, "SKIPPED INITIALISATION gamma_reg=", 0.0
            end if
#endif
            mdl%cost_reg = 0.0_rp
        end if
    else
        mdl%cost_reg = mdl%gamma_reg * sumreg
    end if
    cost = mdl%cost_obs + mdl%cost_reg
#ifndef CPP_ADJ    
    if (mdl%gamma_reg > 1e-12_rp) then
        print *, "COST: Jobs=", mdl%cost_obs
        print *, "      Jreg=", mdl%cost_reg
        print *, " gamma_reg=", mdl%gamma_reg
    end if
#endif
#endif


#ifdef REGUL_THIBAULT
!!! EN TRAVAUX

    !sumreg = ctrl%calc_cost_reg(mdl%msh)!0.0

    !!!  terme de reg => gradient de la bathy-bathy background
    !if (ctrl%Bdemi%n > 0) then
    !    ics = mdl%msh%seg(1)%first_cs
    !    sumreg = ( (ctrl%x(2)-ctrl%x(1)) / mdl%msh%cs(ics+1)%deltademi )**2
    !    do ics = mdl%msh%seg(1)%first_cs+1, mdl%msh%seg(1)%last_cs-1
    !        icsp = ics-mdl%msh%seg(1)%first_cs+2
    !        icsm = ics-mdl%msh%seg(1)%first_cs
    !        sumreg=sumreg+((ctrl%x(icsp)-ctrl%x(icsm))/(2*mdl%msh%cs(ics)%deltademi))**2     
    !    end do
    !    ics = mdl%msh%seg(1)%last_cs
    !    icsp = ics-mdl%msh%seg(1)%first_cs+1
    !    icsm = ics-mdl%msh%seg(1)%first_cs
    !    sumreg=sumreg+((ctrl%x(icsp)-ctrl%x(icsm))/(mdl%msh%cs(ics)%deltademi))**2 
    !else
    !    ics = mdl%msh%seg(1)%first_cs
    !    sumreg = ( ((ctrl%x(2)-ctrl%x0(2))-(ctrl%x(1)-ctrl%x0(1))) / mdl%msh%cs(ics+1)%deltademi )**2
    !    do ics = mdl%msh%seg(1)%first_cs+1, mdl%msh%seg(1)%last_cs-1
    !	     icsp = ics-mdl%msh%seg(1)%first_cs+2
    !        icsm = ics-mdl%msh%seg(1)%first_cs
    !        sumreg=sumreg+((ctrl%x(icsp)-ctrl%x0(icsp)-ctrl%x(icsm)+ctrl%x0(icsm))/(2*mdl%msh%cs(ics)%deltademi))**2  
    !    end do
    !    ics = mdl%msh%seg(1)%last_cs
    !    icsp = ics-mdl%msh%seg(1)%first_cs+1
    !    icsm = ics-mdl%msh%seg(1)%first_cs
    !    sumreg=sumreg+((ctrl%x(icsp)-ctrl%x0(icsp)-ctrl%x(icsm)+ctrl%x0(icsm))/(mdl%msh%cs(ics)%deltademi))**2 
    !end if
 
    !!!  terme de reg => gradient de la bathy
    !ics = mdl%msh%seg(1)%first_cs
    !sumreg = ( (mdl%msh%cs(ics+1)%bathy-mdl%msh%cs(ics)%bathy) / mdl%msh%cs(ics+1)%deltademi )**2
    !do ics = mdl%msh%seg(1)%first_cs+1, mdl%msh%seg(1)%last_cs-1
    !    sumreg = sumreg + ( (mdl%msh%cs(ics+1)%bathy-mdl%msh%cs(ics-1)%bathy) / (2*mdl%msh%cs(ics)%deltademi) )**2     
    !end do
    !ics = mdl%msh%seg(1)%last_cs
    !sumreg = sumreg + ( (mdl%msh%cs(ics)%bathy-mdl%msh%cs(ics-1)%bathy) / (mdl%msh%cs(ics)%deltademi) )**2 

!!! terme de reg => ecart au background
    !if (ctrl%Bdemi%n > 0) then
    !    sumreg = sum((ctrl%x(:))**2)  
    !else
    !    sumreg = sum((ctrl%x(:) - ctrl%x0(:))**2)
    !end if 

    !if (sumreg < 0) then 
    !    write(6,*) "HOLY SHIT, regularization term is negative"
    !end if
    !if (cost < 0) then 
    !    write(6,*) "HOLY SHIT, observation term is negative"
    !end if
    !write(6,*) "HOLY SHIT, the reg is :",sumreg
    !cost = cost + ctrl%coef_reg*sumreg
#endif


#ifdef CPP_ADJ    
    ! Fake operation for Tapenade to ensure that call to cost is made
    cost = sqrt( cost**2 )
#endif
!     cost = sum((obs%est(1, :) - obs%obs(1, :))**2)
    
end subroutine calc_cost
