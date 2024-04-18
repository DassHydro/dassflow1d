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
!> @file update_probes.f90
!> @brief TODO

! #define DEBUG_L1
!> \brief TODO
subroutine update_probes(mdl)
    use m_common
    use m_sw_mono, only: Model
    use m_numeric
#ifdef DEBUG_L1
    use, intrinsic :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param mdl Model
    TYPE(Model), intent(inout) :: mdl
    !=LOCAL VARIABLES==================================================================================================!
    ! Iterator
    integer(ip) :: i
    ! Index of estimation data
    integer(ip) :: idata
    ! Index of observation time
    integer(ip) :: it
    ! Index of probe
    integer(ip) :: iprobe
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of variable
    integer(ip) :: ivar
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
    
    ! Return if there is no probe 
    if (.not. allocated(mdl%probes)) return
    
    ! Loop on probes
    do iprobe = 1, size(mdl%probes)
    
        it = mdl%probes(iprobe)%it
        if (it > size(mdl%probes(iprobe)%t)) cycle
        
        if (mdl%tc >= mdl%probes(iprobe)%t(it)) then
        
            if (it < size(mdl%probes(iprobe)%t)) then
                do while (mdl%tc >= mdl%probes(iprobe)%t(it+1))
                    mdl%probes(iprobe)%it = mdl%probes(iprobe)%it + 1
                    it = it + 1
                    if (it == size(mdl%probes(iprobe)%t)) exit
                end do
            end if
            
            ics = mdl%probes(iprobe)%ics
            
            do ivar = 1, size(mdl%probes(iprobe)%variables)
            
                if (mdl%probes(iprobe)%variables(ivar) == 'Q') then
                
                    mdl%probes(iprobe)%data(ivar, it) = mdl%dof%q(ics)
                
                end if
            
                if (mdl%probes(iprobe)%variables(ivar) == 'H' .or. &
                    mdl%probes(iprobe)%variables(ivar) == 'Z') then
                
                    mdl%probes(iprobe)%data(ivar, it) = mdl%dof%h(ics) + mdl%msh%cs(ics)%bathy
                
                end if
                
            end do

            mdl%probes(iprobe)%it = mdl%probes(iprobe)%it + 1
        
        end if
    
    end do
    
end subroutine update_probes
