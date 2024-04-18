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
!> @file calc_estimations.f90
!> @brief TODO

! #define DEBUG_L1
!> \brief TODO
subroutine calc_estimations(mdl, obs)
    use m_common
    use m_sw_mono, only: Model
    use m_numeric
    use m_obs, only: Observations
#ifdef DEBUG_L1
    use, intrinsic :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param mdl Model
    TYPE(Model), intent(in) :: mdl
    !> @param obs Observations
    type(Observations), intent(inout) :: obs
    !=LOCAL VARIABLES==================================================================================================!
    ! Iterator
    integer(ip) :: i
    ! Index of estimation data
    integer(ip) :: idata
    ! Index of observation time
    integer(ip) :: iobs
    ! Index of station
    integer(ip) :: ista
    ! Index of cross-section
    integer(ip) :: ics
    ! Mean height
    real(rp) :: H_mean
    ! Multiplier
    real(rp) :: mult
    ! Mean discharge
    real(rp) :: q_mean
    ! Width
    real(rp) :: w
    ! Mean width
    real(rp) :: w_mean
    
    ! Loop on stations
    do ista = 1, size(obs%stations)
    
        iobs = obs%stations(ista)%iobs
        if (iobs > size(obs%stations(ista)%t)) cycle
        
!         print *, mdl%tc
!         print *, iobs, size(obs%stations(ista)%t), obs%stations(ista)%t(iobs)
      
        if (mdl%tc >= obs%stations(ista)%t(iobs)) then
        
            if (iobs < size(obs%stations(ista)%t)) then
                do while (mdl%tc >= obs%stations(ista)%t(iobs+1))
!                     print *, "skip", iobs, obs%stations(ista)%t(iobs)
                    obs%stations(ista)%iobs = obs%stations(ista)%iobs + 1
                    iobs = iobs + 1
                    if (iobs == size(obs%stations(ista)%t)) exit
                end do
!                 read(*,*)
            end if
!             print *, "passed"
        
            H_mean = 0.0_rp
            W_mean = 0.0_rp
#ifndef CPP_ADJ
            if (mdl%discharge_estimation) q_mean = 0.0_rp
#endif
            do i = 1, size(obs%stations(ista)%ics)
            
                ics = obs%stations(ista)%ics(i)
!                 print *, ista, i, ics
                H_mean = H_mean + mdl%dof%h(ics) + mdl%msh%cs(ics)%bathy
#ifdef DEBUG_L1                
                if (isnan(H_mean)) then
                    print *, ista, i, ics
                    print *, "H_mean* is nan:", mdl%dof%h(ics), mdl%msh%cs(ics)%bathy
                    read(*,*)
                end if
#endif
                call width(mdl%msh%cs(ics), mdl%dof%h(ics), w)
                w_mean = w_mean + w
#ifndef CPP_ADJ
                if (mdl%discharge_estimation) q_mean = q_mean + mdl%dof%q(ics)
#endif
                
            end do
            mult = 1.0_rp / size(obs%stations(ista)%ics)
#ifdef DEBUG_L1                
            if (isnan(mult) .or. .not. IEEE_IS_FINITE(mult)) then
                print *, "mult is not finite:", ista, size(obs%stations(ista)%ics)
                read(*,*)
            end if
#endif
            H_mean = H_mean * mult
            w_mean = w_mean * mult
            idata = obs%stations(ista)%offset + iobs
!             print *, "Hmean=", H_mean, obs%obs(1, idata)
!             print *, "bathy=", mdl%msh%cs(ics)%bathy
            obs%est(1, idata) = H_mean
            obs%est(2, idata) = w_mean
#ifdef DEBUG_L1                
            if (isnan(H_mean)) then
                print *, "H_mean is nan:", ista, idata, mult
                read(*,*)
            end if
            if (isnan(W_mean)) then
                print *, "W_mean is nan:", ista, idata
                read(*,*)
            end if
#endif
#ifndef CPP_ADJ
            if (mdl%discharge_estimation) then
                q_mean = q_mean * mult
                obs%est(3, idata) = q_mean
            end if
#endif
!             print *, ista, iobs, obs%stations(ista)%offset, idata, H_mean
            
            obs%stations(ista)%iobs = obs%stations(ista)%iobs + 1
        
        end if
    
    end do
    
end subroutine calc_estimations
