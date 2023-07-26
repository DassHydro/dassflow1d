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
!> @file preissmann_time_step.f90
!> @brief This file includes the subroutine for solving the Shallow-Water Equations using Preissmann scheme with LPI 
!> and double sweep methods.


!> \brief Solve the Shallow-Water Equations using Preissmann scheme with LPI and double sweep methods.
subroutine time_loop(mdl, obs)
    use m_common
    use m_mesh
!     use m_sw_mono, only: ImplicitMatrix
    use m_sw_mono, only: Model
    use m_sw_mono, only: reset_probes
!     use m_sw_mono, only: Unknowns
    use m_numeric
    use m_obs
#ifndef CPP_ADJ
    use m_stdout
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param mdl Model
    TYPE(Model), intent(inout) :: mdl
!     !> @param msh Mesh of the model
!     TYPE(Mesh), intent(inout) :: msh
!     !> @param imp Implicit matrix
!     TYPE(ImplicitMatrix), intent(inout) :: imp
!     !> @param dof Unknowns of the model
!     TYPE(Unknowns), intent(inout) :: dof
#ifndef CPP_ADJ
    !> @param obs Observations
    type(Observations), intent(inout), optional :: obs
#else
    !> @param obs Observations
    type(Observations), intent(inout) :: obs
#endif
    !=LOCAL VARIABLES==================================================================================================!
    ! Iterator
    integer(ip) :: i
    ! Index of estimation data
    integer(ip) :: idata
    ! Index of observation time
    integer(ip) :: iobs
    ! Index of result
    integer(ip) :: iout
    ! Index of station
    integer(ip) :: ista
    ! Index of cross-section
    integer(ip) :: ics
    ! Number of results
    integer(ip) :: nout
    ! Number of computed timesteps
    integer(ip) :: nt
    ! Inverse of timestep
    real(rp) :: ddt
    ! Mean height
    real(rp) :: H_mean
    ! Multiplier
    real(rp) :: mult
    ! Width
    real(rp) :: w
    ! Mean width
    real(rp) :: w_mean
#ifndef CPP_ADJ
    ! Temporary string
    character(len=32) :: tmp
#endif
    
    mdl%status = 0
    
#ifndef CPP_ADJ
    ! CHECK-UP
    if (mdl%scheme == "undefined") call abort_solver("Scheme is undefined")
#endif

#ifndef CPP_ADJ
    ! Reset counters
    mdl%warning_counters(:) = 0
    mdl%internal_counters(:) = 0
#endif

    ! Reset observations counters
#ifndef CPP_ADJ
    if (present(obs)) then
        call reset_observations_counters(obs)
    end if
#else
    call reset_observations_counters(obs)
#endif

#ifndef CPP_ADJ
    ! Reset probes
    if (allocated(mdl%probes)) then
        call reset_probes(mdl) 
    end if
#endif

#ifndef CPP_ADJ
    ! Initialisation of results
    if (mdl%dtout > 0.0) then
        nout = (mdl%te - mdl%ts) / mdl%dtout
        if (allocated(mdl%res%t)) then
            if (size(mdl%res%t) /= nout + 1) then
                deallocate(mdl%res%t)
                deallocate(mdl%res%q)
                deallocate(mdl%res%h)
                deallocate(mdl%res%a)
                allocate(mdl%res%t(0:nout))
                allocate(mdl%res%q(mdl%msh%ncs, 0:nout))
                allocate(mdl%res%h(mdl%msh%ncs, 0:nout))
                allocate(mdl%res%a(mdl%msh%ncs, 0:nout))
            end if
        else
            allocate(mdl%res%t(0:nout))
            allocate(mdl%res%q(mdl%msh%ncs, 0:nout))
            allocate(mdl%res%h(mdl%msh%ncs, 0:nout))
            allocate(mdl%res%a(mdl%msh%ncs, 0:nout))
        end if
        mdl%res%t(:) = -1e+99_rp
        mdl%res%q(:, :) = -1e+99_rp
        mdl%res%h(:, :) = -1e+99_rp
        mdl%res%a(:, :) = -1e+99_rp
    end if
#endif

    ! Initialisation
    mdl%tc = mdl%ts
    if (.not. mdl%external_initialisation) then
        call standard_step(mdl, 1.d-4, 1000)
        if (mdl%dtout > 0.0) then
            mdl%res%t(0) = mdl%tc
            mdl%res%h(:, 0) = mdl%dof%h(:)
            mdl%res%q(:, 0) = mdl%dof%q(:)
            mdl%res%a(:, 0) = mdl%dof%a(:)
        end if
        if (mdl%status /= 0) return
    end if
#ifndef CPP_ADJ
    if (mdl%dtout > 0.0) then
        mdl%res%t(0) = mdl%tc
        mdl%res%h(:, 0) = mdl%dof%h(:)
        mdl%res%q(:, 0) = mdl%dof%q(:)
        mdl%res%a(:, 0) = mdl%dof%a(:)
            
        if (len_trim(mdl%output_file) > 0) then
            call init_results(mdl)
            call write_results(mdl, 0)
        end if
        
    end if
#endif
        
    ! Compute estimations
#ifndef CPP_ADJ
    if (present(obs)) then
#endif
        call calc_estimations(mdl, obs)
#ifndef CPP_ADJ
    end if
#endif

    ! Update probes
#ifndef CPP_ADJ
    call update_probes(mdl)
#endif

    ! Update internal counters
#ifndef CPP_ADJ
    call update_internal_counters(mdl)
#endif
   
    ! Time loop
    nt = 0
    do while(mdl%tc < mdl%te)
   
        ! Advance time
        nt = nt + 1
        mdl%tc = mdl%tc + mdl%dt
#ifndef CPP_ADJ
        if (int((mdl%tc - mdl%dt) / mdl%dtout) < int(mdl%tc / mdl%dtout)) then
            if (mdl%print_progress) then
                call print_stdout('dt', (/nt/), (/mdl%tc, mdl%te, (mdl%tc - mdl%ts) / (mdl%te - mdl%ts), mdl%dt/))
            end if
        end if
#endif
        
        ! Solve timestep
        if (mdl%scheme == "preissmann") then
            call preissmann_timestep(mdl, mdl%msh, mdl%imp, mdl%dof, mdl%status)
        end if
	if (mdl%scheme == "implicit_diffusive_wave") then
            call implicit_diffusive_wave(mdl, mdl%msh, mdl%large_grid, mdl%imp, mdl%dof, mdl%status)
        end if

        if (mdl%status /= 0 .and. (.not. mdl%disable_stdout)) then
            print *, "exit:status=", mdl%status
            return
        end if
        
        ! Compute estimations
#ifndef CPP_ADJ
        if (present(obs)) then
#endif
            call calc_estimations(mdl, obs)
#ifndef CPP_ADJ
        end if
#endif
        
        ! Update probes
#ifndef CPP_ADJ
        call update_probes(mdl)
#endif
        
        ! Update internal counters
#ifndef CPP_ADJ
        call update_internal_counters(mdl)
#endif

        ! Store results
#ifndef CPP_ADJ
        if (mdl%dtout > 0.0) then
        
            if (int((mdl%tc - mdl%dt) / mdl%dtout) < int(mdl%tc / mdl%dtout)) then
                iout = int((mdl%tc - mdl%ts) / mdl%dtout)
                mdl%res%t(iout) = mdl%tc
                mdl%res%h(:, iout) = mdl%dof%h(:)
                mdl%res%q(:, iout) = mdl%dof%q(:)
                mdl%res%a(:, iout) = mdl%dof%a(:)
            
                if (len_trim(mdl%output_file) > 0) then
                    call write_results(mdl, iout)
                end if
                
            end if
            
        end if
#endif

    end do
    
#ifndef CPP_ADJ
    ! Display warnings
    if (mdl%warning_counters(1) > 0 .and. (.not. mdl%disable_stdout)) then
        write(tmp, '(I32)') mdl%warning_counters(1)
        print '(3A)', "[ WARNING ] heps correction applied ", trim(adjustl(tmp)), " times"
    end if
#endif
   
end subroutine time_loop

