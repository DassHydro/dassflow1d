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
subroutine calc_cost_and_gradients(cost, grad, mdl, ctrl, obs)
    use m_common
    use m_control
    use m_control_back
    use m_sw_mono, only: Model
    use m_sw_mono_back, only: Model_back
!     use m_numeric
    use m_obs
    use m_obs_back
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param mdl Model
    type(Model), intent(inout) :: mdl
    !> @param ctrl Control
    type(Control), intent(inout) :: ctrl
    !> @param obs Observations
    type(Observations), intent(inout) :: obs
    !> @param dof Unknowns of the model
    real(rp), intent(out) :: cost
    !> @param dof Unknowns of the model
    real(rp), dimension(size(ctrl%x)), intent(out) :: grad
    !=LOCAL VARIABLES==================================================================================================!
    logical :: deallocate_strickler_fields_at_end
    ! Adjoint of control
    type(Control_back) :: ctrl_back
    ! Adjoint of cost
    real(rp) :: cost_back
    ! Adjoint of unknowns
    type(Unknowns_back) :: dof_back
    ! Temporary variable for storing dtout
    real(rp) :: dtout_tmp
    ! Adjoint of Observations
    type(Observations_back) :: obs_back
    ! Adjoint of model
    type(Model_back) :: mdl_back
!     ! Iterator
!     integer(ip) :: i
!     ! Index of estimation data
!     integer(ip) :: idata
!     ! Index of observation time
!     integer(ip) :: iobs
!     ! Index of station
!     integer(ip) :: ista
!     ! Index of cross-section
!     integer(ip) :: ics
!     ! Inverse of timestep
!     real(rp) :: ddt
!     ! Mean height
!     real(rp) :: H_mean
!     ! Multiplier
!     real(rp) :: mult
!     ! Width
!     real(rp) :: w
!     ! Mean width
!     real(rp) :: w_mean
    
#ifndef CPP_ADJ
    ! CHECK-UP
    if (mdl%scheme == "undefined") call abort_solver("Scheme is undefined")
#endif
    
    ! Initialise obs/est arrays if necessary
    if (.not.allocated(obs%obs)) then
        if (.not.allocated(obs%stations)) then
            call abort_solver("No observations defined")
        end if
        call data_from_stations(obs, mdl%discharge_estimation)
    end if
    
    ! Initialise back variables
    call init_back(cost_back, mdl, mdl_back, ctrl, ctrl_back, obs, obs_back)
!     
!     ! Ensure mdl%msh%strickler_fields has a size (required by adjoint)
!     if (allocated(mdl%msh%strickler_fields)) then
!         deallocate_strickler_fields_at_end = .false.
!     else
!         allocate(mdl%msh%strickler_fields(1))
!         deallocate_strickler_fields_at_end = .true.
!     end if
    
    ! Call calc_cost_back
    dtout_tmp = mdl%dtout
    mdl%dtout = -1.0_rp
    call calc_cost_back(cost, cost_back, mdl, mdl_back, ctrl, ctrl_back, obs, obs_back)
    mdl%dtout = dtout_tmp
!     
!     if (deallocate_strickler_fields_at_end) then
!         deallocate(mdl%msh%strickler_fields)
!     end if
    
    ! Transfer gradients
    grad(:) = ctrl_back%x(:)
    
    ! Display warnings
    if (mdl%status /= 0 .and. (.not. mdl%disable_stdout)) then
        print '(A)', "[ WARNING ] simulation failed"
    end if
   
end subroutine calc_cost_and_gradients
