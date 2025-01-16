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

#ifndef CPP_ADJ
! #define DEBUG_L1
! #define DEBUG_L3
! #define DEBUG_DEBORD
! #define OLD_VERSION
#endif
!> \brief Solve the Shallow-Water Equations using Preissmann scheme with LPI and double sweep methods.
subroutine preissmann_timestep(mdl, msh, imp, dof, status)
    use m_common
    use m_mesh
    use m_sw_mono, only: ImplicitMatrix
    use m_sw_mono, only: Model
    use m_sw_mono, only: Unknowns
    use m_numeric
#ifndef CPP_ADJ
    use, intrinsic :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param msh Model
    TYPE(Model), intent(inout) :: mdl
    !> @param msh Mesh of the model
    TYPE(Mesh), intent(inout) :: msh
    !> @param imp Implicit matrix
    TYPE(ImplicitMatrix), intent(inout) :: imp
    !> @param dof Unknowns of the model
    TYPE(Unknowns), intent(inout) :: dof
    !> @param status Status code
    integer(ip) :: status
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of downstream segment
    integer(ip) :: ds_seg
    ! Iterator
    integer(ip) :: i
    ! Index of boundary condition
    integer(ip) :: ibc
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of edge
    integer(ip) :: ie
    ! Index of segment
    integer(ip) :: iseg
    ! Iterator on upstream segments
    integer(ip) :: ius_seg
    ! Offset in non zero entries
    integer(ip) :: off
    ! Index of row
    integer(ip) :: row
    ! Index of upstream segment
    integer(ip) :: us_seg
    ! Some state variables for left cross-section
    real(rp) :: aL, debL, FrL, kL, pL, qL, rhL, uL, zL, wL, ddaL, ddebdzL, dpdzL, dkdzL
    ! Some state variables for right cross-section
    real(rp) :: aR, debR, FrR, kR, pR, qR, rhR, uR, zR, wR, ddaR, ddebdzR, dpdzR, dkdzR
    ! Convection coefficients for left and right cross-sections
    real(rp) :: betaL, betaR
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
    ! Coumpound channel state variables for left cross-section
    real(rp), dimension(3) :: acL, kcL, pcL, rhcL, wcL, dpcdzL
    ! Coumpound channel state variables for right cross-section
    real(rp), dimension(3) :: acR, kcR, pcR, rhcR, wcR, dpcdzR
#endif
#ifdef DEBORD_FORMULA
    ! A0 and A coefficient for Debord Formula
    real(rp) :: A0d, Ad, AdL, AdR
    ! Ratio of hydraulic radiuses (floodplain over channel)
    real(rp) :: r
    ! Ratio of discharges (floodplain over channel)
    real(rp) :: eta
    ! Generalized conveyances
    real(rp), dimension(2) :: debg
#endif
    ! Inverse of timestep
    real(rp) :: ddt
    ! Inverse of gravity acceleration
    real(rp) :: ddg
    ! Inverse of spacing between left and right cross-section
    real(rp) :: ddx
    ! Multiplier
    real(rp) :: mult
    ! Implicit coefficient
    real(rp) :: theta
    ! LPI coefficient
    real(rp) :: thetaLPI
#ifdef DEBUG_L1
    real(rp), dimension(size(dof%h)) :: hprec
    real(rp), dimension(size(dof%h)) :: qprec
    logical :: has_too_big_depth
#endif
    
    intrinsic isnan
    
    real(rp), dimension(msh%ncs) :: DQ
    real(rp), dimension(msh%ncs) :: DZ
    
    
!     real(rp), dimension(4, nseg) :: TA
!     real(rp), dimension(2, nseg) :: TB
!     real(rp), dimension(4, nseg) :: TC
!     real(rp), dimension(2, nseg) :: TD

    real(rp) :: sigma          ! Value of sigma in the preissmann scheme

   ! Definition C coefficients
   real(rp) :: Cb1,Cb2,Cb3          ! Coefficient Cb1,Cb2,Cb3
   real(rp) :: Cc1,Cc4,Cc5          ! Coefficient Cc1,Cc4,Cc5
   real(rp) :: Cd1,Cd4,Cd5          ! Coefficient Cd1,Cd4,Cd5
   real(rp) :: Ce2,Ce3              ! Coefficient Ce2,Ce3
   real(rp) :: Clat1,Clat2          ! Coefficient Clat1,Clat2 
   real(rp) :: Cl1,Cl2,Cl3,Cl4,Cl5  ! Coefficient Cl1,Cl2,Cl3,Cl4,Cl5
   real(rp) :: Cf1,Cf2,Cf3,Cf4,Cf5  ! Coefficient Cf1,Cf2,Cf3,Cf4,Cf5
   real(rp) :: Cg,Ch,Ci,Cj,Ck       ! Coefficient Cg,Ch,Ci,Cj,Ck
   real(rp) :: Ch1,Ch4,Ch5          ! Coefficient Ch1,Ch4,Ch5
   real(rp) :: Cl,Cm,Cn,Co,Cp       ! Coefficient Cl,Cm,Cn,Co,Cp  
   real(rp), dimension(msh%ncs+4) :: cr,cs,ct ! Coefficient Cr, Cs, Ct
   
   ! Definition of G coefficients
   real(rp), dimension(msh%ncs+4) :: Gaj,Gbj,Gcj ! Coefficient Gaj,Gbj,Gcj 
   real(rp), dimension(msh%ncs+4) :: Gdj,Gej,Gfj ! Coefficient Gdj,Gej,Gfj 

   ! Useful computation
   real(rp) :: dx,dpdx             ! Spatial step size and inverse of spatial step size
   
   ! Physical parameters
!    real(rp), dimension(msh%ncs+4) :: q         ! Discharge
!    real(rp), dimension(2) :: qlat      ! Lateral discharge at x and x+1
!    real(rp), dimension(2) :: qlatnext  ! Lateral discharge at x and x+1, at t+1
!    real(rp), dimension(msh%ncs+4) :: z         ! Water elevation
!    real(rp), dimension(msh%ncs+4) :: w         ! Water surface width 
!    real(rp), dimension(msh%ncs+4) :: perimeter ! Wet Perimeter 
!    real(rp), dimension(msh%ncs+4) :: dPdZ      ! dPdZ
!    real(rp), dimension(msh%ncs+4) :: dKdZ      ! dKdZ
!    real(rp), dimension(msh%ncs+4) :: rh        ! Hydraulic radius
!    real(rp), dimension(msh%ncs+4) :: v         ! Velocity flow
!    real(rp), dimension(msh%ncs+4) :: Manning   ! Manning
!    real(rp), dimension(msh%ncs+4) :: debitance ! Debitance
!    real(rp), dimension(msh%ncs+4) :: zp1       ! New Water elevation
  

   integer(ip), dimension(msh%ncs+4) :: mVector              ! Pressure for fluxes computing

   real(rp) :: val0, val,conlim,dcldval

   !Function
   real(rp)  :: dPdZFromH
   real(rp)  :: PerimeterFromH
   real(rp)  :: GetWFromH
   real(rp)  :: HtoS

   real(rp)  :: zaval,qaval,temp,temp1,So

   real(rp)  :: coeff,a,fr,b,coeff1,coeff2,frp
      
   real(rp)  :: flow_ratio

   real(rp)  :: frlpi=0.7_rp
   real(rp)  :: mlpi=10.0_rp
   integer(ip) :: i2, iloc, coefnext
   
   
   integer(ip) :: iic
   real(rp) :: qlat(2)

   !======================================================================================================================!
   !  Update of data
   !======================================================================================================================!
#ifdef DEBUG_L1
    hprec(:) = dof%h(:)
    qprec(:) = dof%q(:)
    do ics = 1, size(msh%cs)
        if (isnan(dof%h(ics))) then
            print *, "FIRST CHECK;h is NaN:", ics
            read(*,*)
        end if
    end do
#endif

   call update_all_levels(msh, dof%h)
   
#ifdef DEBUG_L1
    do ics = 1, size(msh%cs)
        if (isnan(dof%h(ics))) then
            print *, "SECOND CHECK;h is NaN:", ics
            read(*,*)
        end if
    end do
#endif
            
    ! Compute inverse of timestep and gravity
    ddt = 1.0_rp / mdl%dt
    ddg = 1.0_rp / mdl%gravity
    
    ! Compute transfer coefficients (GA, GB, GC, GD, GE, GF)
    iic = 1
    do iseg = 1, msh%nseg
    
!         print *, "BIEF", iseg
        
        ! Retrieve state variables for first left cross-section
        ie = msh%seg(iseg)%first_cs
        zL = dof%h(ie) + msh%cs(ie)%bathy
        qL = dof%Q(ie)
        if (msh%strickler_type_code == strickler_type_average_conveyance) then
#ifdef AVERAGE_CONVEYANCE
            call areas_compound_channel(msh%cs(ie), dof%h(ie), acL)
            call perimeters_compound_channel(msh%cs(ie), dof%h(ie), pcL)
            call widths_compound_channel(msh%cs(ie), dof%h(ie), wcL)
            kcL(:) = msh%cs(ie)%strickler_params(:)
            call dpdz_compound_channel(msh%cs(ie), dof%h(ie), dpcdzL)
            do i = 1, 3
                dpcdzL(i) = min(100.0_rp, dpcdzL(i))
                if (pcL(i) > 1.e-12_rp) then
                    rhcL(i) = acL(i) / pcL(i)
                else
                    rhcL(i) = 0.d0
                end if
            end do
            debL = kcL(1) * acL(1) * rhcL(1)**d2p3 + &
                   kcL(2) * acL(2) * rhcL(2)**d2p3 + &
                   kcL(3) * acL(3) * rhcL(3)**d2p3
            ddebdzL = d1p3 * (kcL(1) * rhcL(1)**d2p3 * (5.0_rp * wcL(1) - 2.0_rp * rhcL(1) * dpcdzL(1)) + &
                              kcL(2) * rhcL(2)**d2p3 * (5.0_rp * wcL(2) - 2.0_rp * rhcL(2) * dpcdzL(2)) + &
                              kcL(3) * rhcL(3)**d2p3 * (5.0_rp * wcL(3) - 2.0_rp * rhcL(3) * dpcdzL(3)))
            aL = sum(acL)
            pL = sum(pcL)
            wL = sum(wcL)
            dkdzL = 0.0_rp
            betaL = 1.0_rp
#endif
#ifdef DEBORD_FORMULA
        else if (msh%strickler_type_code == strickler_type_debord) then

            ! Retrieve compound state variables
            kcL(:) = msh%cs(ie)%strickler_params(:)
!             print *, "LEFT:", ie, loc(msh%cs(ie)), msh%cs(ie)%level_heights
            call areas_compound_channel(msh%cs(ie), dof%h(ie), acL)
            call perimeters_compound_channel(msh%cs(ie), dof%h(ie), pcL)
            call widths_compound_channel(msh%cs(ie), dof%h(ie), wcL)
            call dpdz_compound_channel(msh%cs(ie), dof%h(ie), dpcdzL)
            
            ! Compute variables for floodplain (composition of right and left overbanks)
            if (pcL(1) + pcL(3) > zerom) then
                kcL(1) = (kcL(1) * pcL(1) + kcL(3) * pcL(3)) / (pcL(1) + pcL(3))
            else
                kcL(1) = 0.5_rp * (kcL(1) + kcL(3))
            endif
            pcL(1) = pcL(1) + pcL(3)
            acL(1) = acL(1) + acL(3)
            wcL(1) = wcL(1) + wcL(3)
            dpcdzL(1) = dpcdzL(1) + dpcdzL(3)

            ! Compute hydraulic radiuses and limit dPdZ
            do i = 1, 2
                dpcdzL(i) = min(mdl%dpdz_cutoff, dpcdzL(i))
                if (pcL(i) > tinym) then
                    rhcL(i) = acL(i) / pcL(i)
                else
                    rhcL(i) = 0.d0
                end if
            end do
            
            ! Compute A0 and A Debord
            A0d = 0.9_rp * (kcL(1) / kcL(2))**d1p6
            if (rhcL(1) > tinym) then
                r = rhcL(1) / rhcL(2)
                if (r > 0.3_rp) then
                    Ad = A0d
                else
                    Ad = demi * (1.0_rp - A0d) * cos(pi * r / 0.3_rp) + demi * (1.0_rp + A0d)
                end if
            else
                Ad = 1.0_rp
            end if
            AdL = Ad
            
            ! Update Stricklers before conveyance computation
            if (acL(2) > tinym) then
                if (1.0_rp + acL(1) / acL(2) * (1.0_rp - Ad**2) > tinym) then
                    kcL(1) = sqrt(1.0_rp + acL(1) / acL(2) * (1.0_rp - Ad**2)) * kcL(1)
                    kcL(2) = Ad * kcL(2)
                else
! #ifndef CPP_ADJ
!                 print *, "[ SIMULATION FAILED ] Debord operand is negative"
! #endif
                open(199, file="error.txt")
                write(199, '(2(A,I6))') "Negative Debord operand detected on segment ", iseg-1, ", cross-section ", ics
                close(199)

                status = 22
                return
            end if
                
            end if
            
            ! Compute generalized conveyances
            debg(1) = kcL(1) * acL(1) * rhcL(1)**d2p3
            debg(2) = kcL(2) * acL(2) * rhcL(2)**d2p3
                
            ! Compute total conveyance
            debL = debg(1) + debg(2)
            ddebdzL = d1p3 * (kcL(1) * rhcL(1)**d2p3 * (5.0_rp * wcL(1) - 2.0_rp * rhcL(1) * dpcdzL(1)) + &
                              kcL(2) * rhcL(2)**d2p3 * (5.0_rp * wcL(2) - 2.0_rp * rhcL(2) * dpcdzL(2)))
!             print *, "debL=", debL
!             print *, "++ kminL*=", kcL(2), ", kmajL*=", kcL(1)
!             print *, "++ aminL=", acL(2), ", amajL=", acL(1)
!             print *, "++ rhminL=", rhcL(2), ", rhmajL=", rhcL(1)

            ! Compute beta for the convection term
            aL = acL(1) + acL(2)
            betaL = 1.0_rp
            if( acL(1) > tinym ) then
                eta = debg(2) / debg(1)
                betaL = (eta**2 / acL(2) + 1.0_rp / acL(1)) * aL / (1.0_rp + eta)**2
            end if
            
            ! Compute other total variables
            pL = pcL(1) + pcL(2)
            wL = wcL(1) + wcL(2)
            rhL = aL / pL
            dkdzL = 0.0_rp
#endif
        else
            call area(msh%cs(ie), dof%h(ie), aL)
            call perimeter(msh%cs(ie), dof%h(ie), pL)
            call width(msh%cs(ie), dof%h(ie), wL)
            call strickler(msh%cs(ie), dof%h(ie), kL)
            call dpdz(msh%cs(ie), dof%h(ie), dpdzL)
            dpdzL = min(100.0_rp, dpdzL)
            call dkdz(msh%cs(ie), dof%h(ie), dkdzL)
            rhL = aL / pL
            debL = kL * aL * rhL**d2p3
            ddebdzL = dkdzL * aL * rhL**d2p3 + &
                      d1p3 * (kL * rhL**d2p3) * (5.0_rp * wL - 2.0_rp * rhL * dpdzL)
            betaL = 1.0_rp
        end if
!         betaL = 1.0_rp
        uL = qL / aL
        frL = sqrt(dof%q(ie)**2 * wL / (mdl%gravity * aL**3))
        ddaL = 1.0_rp / aL
        
        ! Loop on edges of segments
        row = imp%seg_offsets(iseg)
        do ie = msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs - 1
        
            row = row + 1
            
            ! Retrieve state variables for right cross-section
            zR = dof%h(ie+1) + msh%cs(ie+1)%bathy
            qR = dof%Q(ie+1)
            if (msh%strickler_type_code == strickler_type_average_conveyance) then
#ifdef AVERAGE_CONVEYANCE
                call areas_compound_channel(msh%cs(ie+1), dof%h(ie+1), acR)
                call perimeters_compound_channel(msh%cs(ie+1), dof%h(ie+1), pcR)
                call widths_compound_channel(msh%cs(ie+1), dof%h(ie+1), wcR)
                kcR(:) = msh%cs(ie+1)%strickler_params(:)
                call dpdz_compound_channel(msh%cs(ie+1), dof%h(ie+1), dpcdzR)
                do i = 1, 3
                    dpcdzR(i) = min(100.0_rp, dpcdzR(i))
                    if (pcR(i) > 1.e-12_rp) then
                        rhcR(i) = acR(i) / pcR(i)
                    else
                        rhcR(i) = 0.d0
                    end if
                end do
                debR = kcR(1) * acR(1) * rhcR(1)**d2p3 + &
                       kcR(2) * acR(2) * rhcR(2)**d2p3 + &
                       kcR(3) * acR(3) * rhcR(3)**d2p3
                ddebdzR = d1p3 * (kcR(1) * rhcR(1)**d2p3 * (5.0_rp * wcR(1) - 2.0_rp * rhcR(1) * dpcdzR(1)) + &
                                  kcR(2) * rhcR(2)**d2p3 * (5.0_rp * wcR(2) - 2.0_rp * rhcR(2) * dpcdzR(2)) + &
                                  kcR(3) * rhcR(3)**d2p3 * (5.0_rp * wcR(3) - 2.0_rp * rhcR(3) * dpcdzR(3)))
                aR = sum(acR)
                pR = sum(pcR)
                wR = sum(wcR)
                dkdzR = 0.0_rp
                betaR = 1.0_rp
#endif
#ifdef DEBORD_FORMULA
        else if (msh%strickler_type_code == strickler_type_debord) then

            ! Retrieve compound state variables
            kcR(:) = msh%cs(ie+1)%strickler_params(:)
!             print *, "RIGHT:", ie+1, loc(msh%cs(ie+1)), msh%cs(ie+1)%level_heights
            call areas_compound_channel(msh%cs(ie+1), dof%h(ie+1), acR)
            call perimeters_compound_channel(msh%cs(ie+1), dof%h(ie+1), pcR)
            call widths_compound_channel(msh%cs(ie+1), dof%h(ie+1), wcR)
            call dpdz_compound_channel(msh%cs(ie+1), dof%h(ie+1), dpcdzR)
            
            ! Compute variables for floodplain (composition of right and left overbanks)
            if (pcR(1) + pcR(3) > zerom) then
                kcR(1) = (kcR(1) * pcR(1) + kcR(3) * pcR(3)) / (pcR(1) + pcR(3))
            else
                kcR(1) = 0.5_rp * (kcR(1) + kcR(3))
            endif
            pcR(1) = pcR(1) + pcR(3)
            acR(1) = acR(1) + acR(3)
            wcR(1) = wcR(1) + wcR(3)
            dpcdzR(1) = dpcdzR(1) + dpcdzR(3)

            ! Compute hydraulic radiuses and limit dPdZ
            do i = 1, 2
                dpcdzR(i) = min(mdl%dpdz_cutoff, dpcdzR(i))
                if (pcR(i) > tinym) then
                    rhcR(i) = acR(i) / pcR(i)
                else
                    rhcR(i) = 0.d0
                end if
            end do
            
            ! Compute A0 and A Debord
            A0d = 0.9_rp * (kcR(1) / kcR(2))**d1p6
            if (rhcR(1) > tinym) then
                r = rhcR(1) / rhcR(2)
                if (r > 0.3_rp) then
                    Ad = A0d
                else
                    Ad = demi * (1.0_rp - A0d) * cos(pi * r / 0.3_rp) + demi * (1.0_rp + A0d)
                end if
            else
                Ad = 1.0_rp
            end if
            AdR = Ad
            
            ! Update Stricklers before conveyance computation
            if (acR(2) > tinym) then
                if (1.0_rp + acR(1) / acR(2) * (1.0_rp - Ad**2) > tinym) then
                    kcR(1) = sqrt(1.0_rp + acR(1) / acR(2) * (1.0_rp - Ad**2)) * kcR(1)
                    kcR(2) = Ad * kcR(2)
                else
! #ifndef CPP_ADJ
!                     print *, "[ SIMULATION FAILED ] Debord operand is negative"
! #endif
                    open(199, file="error.txt")
                    write(199, '(2(A,I6))') "Negative Debord operand detected on segment ", iseg-1, ", cross-section ", ie+1
                    close(199)
                    status = 22
                    return
                end if
            end if
            
            ! Compute generalized conveyances
            debg(1) = kcR(1) * acR(1) * rhcR(1)**d2p3
            debg(2) = kcR(2) * acR(2) * rhcR(2)**d2p3
                
            ! Compute total conveyance
            debR = debg(1) + debg(2)
            ddebdzR = d1p3 * (kcR(1) * rhcR(1)**d2p3 * (5.0_rp * wcR(1) - 2.0_rp * rhcR(1) * dpcdzR(1)) + &
                              kcR(2) * rhcR(2)**d2p3 * (5.0_rp * wcR(2) - 2.0_rp * rhcR(2) * dpcdzR(2)))

            ! Compute beta for the convection term
            aR = acR(1) + acR(2)
            betaR = 1.0_rp
            if( acR(1) > tinym ) then
                eta = debg(2) / debg(1)
                betaR = (eta**2 / acR(2) + 1.0_rp / acR(1)) * aR / (1.0_rp + eta)**2
            end if
            
            ! Compute other total variables
            pR = pcR(1) + pcR(2)
            wR = wcR(1) + wcR(2)
            rhR = aR / pR
            dkdzR = 0.0_rp
#endif
            else
                call area(msh%cs(ie+1), dof%h(ie+1), aR)
                call perimeter(msh%cs(ie+1), dof%h(ie+1), pR)
                call width(msh%cs(ie+1), dof%h(ie+1), wR)
                call strickler(msh%cs(ie+1), dof%h(ie+1), kR)
                call dpdz(msh%cs(ie+1), dof%h(ie+1), dpdzR)
                dpdzR = min(100.0_rp, dpdzR)
                call dkdz(msh%cs(ie+1), dof%h(ie+1), dkdzR)
                rhR = aR / pR
                debR = kR * aR * rhR**d2p3
                ddebdzR = dkdzR * aR * rhR**d2p3 + &
                          d1p3 * (kR * rhR**d2p3) * (5.0_rp * wR - 2.0_rp * rhR * dpdzR)
                betaR = 1.0_rp
            end if
!             betaR = 1.0_rp
            uR = qR / aR
            frR = sqrt(dof%q(ie)**2 * wL / (mdl%gravity * aR**3))
            ddaR = 1.0_rp / aR
#ifdef DEBUG_L1
            if (isNan(debR)) then
                print *, "debR is NaN:", iseg, ie
                print *, "Detail:", kR, aR, rhR, dof%h(ie+1)
            end if
            if (.not. IEEE_IS_FINITE(debR)) then
                print *, "debR is INF:", iseg, ie
                print *, "Detail:", kR, aR, rhR, dof%h(ie+1)
            end if
#endif

#ifndef CPP_ADJ
        ! Check for NaN value(s)
        if (isnan(wL) .or. isnan(aL) .or. isnan(rhL) .or. isnan(qL) .or. &
            isnan(wR) .or. isnan(aR) .or. isnan(rhR)) then
#ifdef DEBUG_L1
            print *, "code 20:", ie
            print *, dof%h(ie), dof%h(ie+1)
            print *, wL, wR
            print *, aL, aR
            print *, rhL, rhR
            read(*,*)
#endif
            status = 20
            open(199, file="error.txt")
            write(199, '(3(A,I6))') "NaN internal variable detected on segment ", iseg-1, ", cross-sections ", ie, &
                                    " and ", ie+1
            write(199, '(A,F6.2,1X,F6.2)') "h:", dof%h(ie), dof%h(ie+1)
            write(199, '(A,F6.2,1X,F6.2)') "b:", msh%cs(ie)%bathy, msh%cs(ie+1)%bathy
            write(199, '(A,F8.2,1X,F8.2)') "w:", wL, wR
            write(199, '(A,E12.5,1X,E12.5)') "a:", aL, aR
            write(199, '(A,F6.2,1X,F6.2)') "r:", rhL, rhR
            write(199, '(A)') ""
            write(199, '(A)') "LEFT:"
            do i = 1, size(msh%cs(ie)%level_heights)
                write(199, '(I4.4,1X,F6.2,1X,F8.2,1X,E12.5)') i, msh%cs(ie)%level_heights(i), msh%cs(ie)%level_widths(i), msh%cs(ie)%area_cum(i)
            end do
            write(199, '(A)') ""
            write(199, '(A)') "RIGHT:"
            do i = 1, size(msh%cs(ie+1)%level_heights)
                write(199, '(I4.4,1X,F6.2,1X,F8.2,1X,E12.5)') i, msh%cs(ie+1)%level_heights(i), msh%cs(ie+1)%level_widths(i), msh%cs(ie+1)%area_cum(i)
            end do
            close(199)
            return
        end if
#endif
!             print *, "EDGE:", ie
!             print *, "K:", kL, kR
!             print *, "P:", pL, pR
!             print *, "dP:", dpdzL, dpdzR
!             print *, "dK:", dkdzL, dkdzR


            ! Compute inverse of spacing and timestep
            dx = msh%cs(ie+1)%deltademi                       !TODO CHECK !!!
            ddx = 1.0_rp / dx
            
            ! Compute CG, CH, CI, CJ, CK (Mass Conservation equation)
            theta = mdl%theta_preissmann
            CG = theta
            CH = 0.25_rp * dx * ddt * (wR + wL)
            CI = CG
            CJ = -CH
            CK = -(qR - qL)
#ifdef DEBUG_L1
            if (isNan(CH)) then
                print *, "CH is NaN:", iseg, ie
                print *, "Detail:", dx, ddt, wR, wL
                print *, "Detail+:", dof%h(ie), dof%h(ie+1)
                read(*,*)
            end if
            if (isNan(CK)) then
                print *, "CK is NaN:", iseg, ie
                print *, "Detail:", qR, qL
            end if
#endif
              
            ! Compute theta and thetaLPI
            if (frL < mdl%frLPI) then
              thetaLPI = 1.0_rp - (FrL / mdl%frLPI)**mdl%mLPI
            else
              thetaLPI = 0.0_rp
            end if
            if (frR < mdl%frLPI) then
              thetaLPI = min(thetaLPI, 1.0_rp - (FrR / mdl%frLPI)**mdl%mLPI)
            else
              thetaLPI = 0.0_rp
            end if
!             print *, iseg, ie, frL, frR
            
            ! Compute CB1, CB2, CB3 (Momentum Conservation equation, term I)
            CB1 = demi * (qR * abs(qR) + qL * abs(qL))
            CB2 = theta * qR
            CB3 = theta * qL
#ifdef DEBUG_L1
            if (isNan(CB1) .or. .not. IEEE_IS_FINITE(CB1)) then
                print *, "CB1 is not finite:", iseg, ie
                print *, "Detail:", qR, qL
            end if
#endif
            
            ! Compute CC1, CC4, CC5 (Momentum Conservation equation, term II)
! #ifdef DEBORD_FORMULA
!             CC1 = demi * (debR**2 + debL**2)
!             CC4 = theta * debR * (d1p3 * kcR * beta_debord * rcR**d2p3 * (5.0_rp * wcR - 2.0_rp * rcR * dpcdzR) + &
!                                   kfR * (0.5_rp * duR / sqrtuR * rfR**d2p3 + &
!                                          sqrtuR * (d2p3 * wfR / (afR**d1p3 * pfR**d2p3 - d2p3 * (rfR / pfR)**d2p3 * dpfdzR))))
!             CC5 = theta * debL * (d1p3 * kcL * beta_debord * rcL**d2p3 * (5.0_rp * wcL - 2.0_rp * rcL * dpcdzL) + &
!                                   kfL * (0.5_rp * duR / sqrtuR * rfL**d2p3 + &
!                                          sqrtuR * (d2p3 * wfL / (afL**d1p3 * pfL**d2p3 - d2p3 * (rfL / pfL)**d2p3 * dpfdzL))))
! #elif defined AVERAGE_CONVEYANCE
!             CC1 = demi * (debR**2 + debL**2)
!             CC4 = theta * debR * (d1p3 * kR(1) * rR(1)**d2p3 * (5.0_rp * wcR(1) - 2.0_rp * rR(1) * dpdzR(1)) + &
!                                   d1p3 * kR(2) * rR(2)**d2p3 * (5.0_rp * wcR(2) - 2.0_rp * rR(2) * dpdzR(2)) + &
!                                   d1p3 * kR(3) * rR(3)**d2p3 * (5.0_rp * wcR(3) - 2.0_rp * rR(3) * dpdzR(3)))
!             CC5 = theta * debL * (d1p3 * kL(1) * rL(1)**d2p3 * (5.0_rp * wcL(1) - 2.0_rp * rL(1) * dpdzL(1)) + &
!                                   d1p3 * kL(2) * rL(2)**d2p3 * (5.0_rp * wcL(2) - 2.0_rp * rL(2) * dpdzL(2)) + &
!                                   d1p3 * kL(3) * rL(3)**d2p3 * (5.0_rp * wcL(3) - 2.0_rp * rL(3) * dpdzL(3)))
! 
! #else
!             CC1 = demi * (debR**2 + debL**2)
!             CC4 = theta * debR * (dkdzR * aR * rhR**d2p3 + &
!                                   d1p3 * (kR * rhR**d2p3) * (5.0_rp * wR - 2.0_rp * rhR * dpdzR))
!             CC5 = theta * debL * (dkdzL * aL * rhL**d2p3 + &
!                                   d1p3 * (kL * rhL**d2p3) * (5.0_rp * wL - 2.0_rp * rhL * dpdzL))
! #endif
            CC1 = demi * (debR**2 + debL**2)
            CC4 = theta * debR * ddebdzR
            CC5 = theta * debL * ddebdzL
#ifdef DEBUG_L1
            if (isNan(CC1)) then
                print *, "CC1 is NaN:", iseg, ie
                print *, "Detail:", debR, debL
            end if
            if (.not. IEEE_IS_FINITE(CC1)) then
                print *, "CC1 is INF:", iseg, ie
                print *, "Detail:", debR, debL
            end if
            if (isNan(CC4)) then
                print *, "CC4 is NaN:", iseg, ie
                print *, "Detail:", debR, dkdzR, aR, rHR, kR, wL, dpdzR
            end if
            if (isNan(CC5)) then
                print *, "CC5 is NaN:", iseg, ie
                print *, "Detail:", debL, dkdzL, aL, rHL, kL, wL, dpdzL
            end if
#endif
            
            ! Compute CD1, CD4, CD5 (Momentum Conservation equation, term III)
            CD1 = ddx * (zR - zL)
            CD4 = theta * ddx
            CD5 = -CD4
#ifdef DEBUG_L1
            if (isNan(CD1)) then
                print *, "CD1 is NaN:", iseg, ie
                print *, "Detail:", ddx, zR, zL
            end if
#endif
            
            ! Compute CH1, CH4, CH5 (Momentum Conservation equation, term IV)
            CH1 = demi * ddg * (ddaR + ddaL)
            CH4 = -demi * theta * wR * ddaR**2 * ddg
            CH5 = -demi * theta * wL * ddaL**2 * ddg
#ifdef DEBUG_L1
            if (isNan(CH1)) then
                print *, "CH1 is NaN:", iseg, ie
                print *, "Detail:", ddg, ddaR, ddaL
            end if
#endif
            
            ! Compute CE2, CE3 (Momentum Conservation equation, term V)
            CE2 = demi * thetaLPI * ddt
            CE3 = CE2
            
            ! Compute CF1, CF2, CF3, CF4, CF5 (Momentum Conservation equation, term VI)
            CF1 = thetaLPI * ddx * (betaR * qR * uR - betaL * qL * uL)
            CF2 = thetaLPI * ddx * two * theta * betaR * uR
            CF3 = -thetaLPI * ddx * two * theta * betaL * uL
            CF4 = -thetaLPI * ddx * theta * betaR * uR**2 * wR
            CF5 = -thetaLPI * ddx * theta * betaL * uL**2 * wL
#ifdef DEBUG_L1
            if (isNan(CF1) .or. .not. IEEE_IS_FINITE(CF1)) then
                print *, "CF1 is not finite:", iseg, ie
                print *, "Detail:", thetaLPI, qR, uR, qL, uL, ddx
            end if
            if (isNan(CF2) .or. .not. IEEE_IS_FINITE(CF2)) then
                print *, "CF2 is not finite:", iseg, ie
                print *, "Detail:", thetaLPI, theta, uR, ddx
            end if
            if (isNan(CF3) .or. .not. IEEE_IS_FINITE(CF3)) then
                print *, "CF3 is not finite:", iseg, ie
                print *, "Detail:", thetaLPI, uL, ddx
            end if
            if (isNan(CF4) .or. .not. IEEE_IS_FINITE(CF4)) then
                print *, "CF4 is not finite:", iseg, ie
                print *, "Detail:", thetaLPI, uR, wR, ddx
            end if
            if (isNan(CF5) .or. .not. IEEE_IS_FINITE(CF5)) then
                print *, "CF5 is not finite:", iseg, ie
                print *, "Detail:", thetaLPI, uL, wL, ddx
            end if
#endif
            
            ! TODO QLAT
            if (allocated(mdl%ic) .and. mdl%nic > 0) then
                if (iic <= size(mdl%ic)) then
                    if (iseg == mdl%ic(iic)%iseg .and. ie == mdl%ic(iic)%ie) then
                        mult = 1.0 / mdl%ic(iic)%ledge
                        qlat(1) = linear_interp(mdl%ic(iic)%ts%t, mdl%ic(iic)%ts%y, mdl%tc-mdl%dt) * mult
                        qlat(2) = linear_interp(mdl%ic(iic)%ts%t, mdl%ic(iic)%ts%y, mdl%tc) * mult
                        CK = CK + dx * (theta * qlat(2) + (1.0 - theta) * qlat(1))
!                         print *, "Ck+=", dx * (theta * qlat(2) + (1.0 - theta) * qlat(1))
!                         print *, "++", dx,  theta, qlat(2), 1.0 - theta, qlat(1)
                        CL1 = thetaLPI * 0.5 * (uL + uR) * (theta * qlat(2) + (1.0 - theta) * qlat(1))
                        mult = thetaLPI * theta * 0.5 * qlat(2)
                        CL2 = mult / aR
                        CL3 = -mult * uR * wR / aR
                        CL4 = mult / aL
                        CL5 = -mult * uL * wL / aL
!                         print *, "Cl:", Cl1, Cl2, Cl3, Cl4, Cl5
!                         print *, "++", thetaLPI, uL, uR
!                         print *, "geoQ_L", qL, aL, dof%h(ie)
!                         print *, "geoQ_R", qR, aR, dof%h(ie+1)
!                         read(*,*)
                        CF1 = CF1 - CL1
                        CF2 = CF2 - CL2
                        CF3 = CF3 - CL3
                        CF4 = CF4 - CL4
                        CF5 = CF5 - CL5
                        iic = iic + 1
                    end if
                end if
            end if
            
! ! !      select case(lat_inflow)
! ! !        case ('1')
! ! ! 
! ! !        if ( ANY( bc%hyd_lat%loc == ie ) ) then
! ! !             iloc=iloc+1
! ! !             qlat(1)      = dof%qlat(1,iloc)                          ! Update of qlat
! ! !             qlatnext(1)  = dof%qlat(2,iloc)                        ! Update of qlat t+1
! ! !             
! ! !             
! ! !             qlat(2)      = 0.0_rp
! ! !             qlatnext(2)  = 0.0_rp
! ! !             coefnext=1
! ! !             if (iloc < size(bc%hyd_lat%loc)) then
! ! !                if ( bc%hyd_lat%loc(iloc+1) == ie+1 ) then
! ! !                   qlat(2)      = dof%qlat(1,iloc+1)                          ! Update of qlat x+1
! ! !                   qlatnext(2)  = dof%qlat(2,iloc+1)                      ! Update of qlat x+1, t+1
! ! !                   coefnext=2
! ! !                endif
! ! !             endif
! ! ! 
! ! ! ! 
! ! ! !       !!!Mass conservation
! ! !         Ck    = Ck + sigma * ( qlatnext(2) + qlatnext(1) )/coefnext + ( 1 - sigma ) * ( qlat(2) + qlat(1) )/coefnext
! ! ! ! 	!!!Momentum conservation
! ! ! ! 	
! ! !  	flow_ratio = 0.1 ! To discretize later : speed_ratio = Qlat/(Q+Qlat)
! ! ! ! 	
! ! ! ! 	!!If ulat is given    
! ! !  	!Clat1 = sigma/2*(ulatnext(ie+1)+qlatnext(ie)) + (1-sigma)/2*(qlat(ie+1)+qlat(ie))
! ! !  	Clat2 = ( sigma / coefnext * ( qlatnext(2) + qlatnext(1) ) + ( 1 - sigma ) / coefnext * ( qlat(2) + qlat(1) ) ) * flow_ratio
! ! !  	!Cp    = -(Cb1+Cc1*(Cd1+Ch1*(Cf1-Clat1*Clat2)))
! ! ! ! 	
! ! ! ! 	!!If ulat is a function of Q and S
! ! !  	Cl1 =  coeff * (q(ie+1) / dof%A(ie+1) - q(ie) / dof%A(ie) ) / 2
! ! !  	Cl2 =  coeff * sigma / 2 / dof%A(ie+1)
! ! !  	Cl3 = -coeff * sigma / 2 / dof%A(ie)
! ! !  	Cl4 = -coeff * sigma / 2 * w(ie+1) * v(ie+1) / dof%A(ie+1)
! ! !  	Cl5 =  coeff * sigma / 2 * w(ie  ) * v(ie  ) / dof%A(ie)
! ! ! ! 
! ! !  	Cf1 = Cf1 - Cl1*Clat2
! ! !  	Cf2 = Cf2 - Cl2*Clat2
! ! !  	Cf3 = Cf3 - Cl3*Clat2
! ! !  	Cf4 = Cf4 - Cl4*Clat2
! ! !  	Cf5 = Cf5 - Cl5*Clat2
! ! ! ! 
! ! ! ! 
! ! !        end if
! ! !      end select  
            
            ! Compute CL, CM, CN, CO, CP (Momentum Conservation equation, all terms)
            CL = CB2 + CC1 * CH1 * (CE2 + CF2)
            CM = CC4 * (CD1 + CH1 * CF1) + CC1 * (CD4 + CH1 * CF4 + CH4 * CF1)
            CN = -(CB3 + CC1 * CH1 * (CE3 + CF3))
            CO = -(CC5 * (CD1 + CH1 * CF1) + CC1 * (CD5 + CH1 * CF5 + CH5 * CF1))
            CP = -(CB1 + CC1 * (CD1 + CH1 * CF1))
#ifdef DEBUG_DEBORD
            if (acL(1) > zerom .or. acR(1) > zerom) then
                print *, "CS(L)", ie
                print *, "-- Q:", dof%Q(ie)
                print *, "-- Z:", msh%cs(ie)%bathy + dof%h(ie)
                print *, "-- h:", dof%h(ie)
                print *, "-- W:", wcL(2), wcL(1), wL
                print *, "CS(R)", ie+1
                print *, "-- Q:", dof%Q(ie+1)
                print *, "-- Z:", msh%cs(ie+1)%bathy + dof%h(ie+1)
                print *, "-- W:", wcR(2), wcR(1), wR
                print *, "-- h:", dof%h(ie+1)
    !             read(*,*)
                print *, "debL=", debL
                print *, "++ kminL*=", kcL(2), ", kmajL*=", kcL(1)
                print *, "++ aminL=", acL(2), ", amajL=", acL(1)
                print *, "++ rhminL=", rhcL(2), ", rhmajL=", rhcL(1)
                print *, "++ AdL=", AdL
                print *, "debR=", debR
                print *, "++ kminR*=", kcR(2), ", kmajR*=", kcR(1)
                print *, "++ aminR=", acR(2), ", amajR=", acR(1)
                print *, "++ rhminR=", rhcR(2), ", rhmajR=", rhcR(1)
                print *, "++ AdR=", AdR
                print *, "CdM:"
                print *, CG, CH, CI, CJ, CK
                print *, "  B:", CB1, CB2, CB3
                print *, "  C:", CC1, CC4, CC5
                print *, "  D:", CD1, CD4, CD5
                print *, "  E:", CE2, CE3
                print *, "  F:", CF1, CF2, CF3, CF4, CF5
                print *, "  H:", CH1, CH4, CH5
                print *, "QdM:"
                print *, CL, CM, CN, CO, CP
                read(*,*)
            end if
#endif
#ifdef DEBUG_L1
            if (isnan(CL) .or. .not. IEEE_IS_FINITE(CL)) then
                print *, "CL is not finite", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                print *, "--details:", CB2, "+", CC1, "*", CH1, "*", "(", CE2, "+", CF2, ")"
                read(*,*)
            end if
            if (isnan(CM) .or. .not. IEEE_IS_FINITE(CM)) then
                print *, "CM is not finite", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                print *, "--details:", CC4, "*(", CD1, "+", CH1, "*", CF1, ")+", CC1, "*(",&
                                       CD4, "+", CH1, "*", CF4, "+", CH4, "*", CF1, ")"
                read(*,*)
            end if
            if (isnan(CN) .or. .not. IEEE_IS_FINITE(CN)) then
                print *, "CN is not finite", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                print *, "--details: -(", CB3, "+", CC1, "*", CH1, "*(", CE3, "+", CF3, "))"
                read(*,*)
            end if
            if (isnan(CO) .or. .not. IEEE_IS_FINITE(CO)) then
                print *, "CO is not finite", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                print *, "--details: -(", CC5, "*(", CD1, "+", CH1, "*", CF1, ")+", CC1, "*(",&
                                       CD5, "+", CH1, "*", CF5, "+", CH5, "*", CF1, "))"
                read(*,*)
            end if
            if (isnan(CP) .or. .not. IEEE_IS_FINITE(CP)) then
                print *, "CP is not finite", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                print *, "--details: -(", CB1, "+", CC1, "*(", CD1, "+", CH1, "*", CF1, "))"
                print *, "Detail(CB1):", qR, qL
                print *, "Detail(Q):", dof%q(:)
                print *, "Detail(h):", dof%h(:)
                !read(*,*)
            end if
#endif
            
            ! Compute transfer coefficients (upstream->downstream)
            mult = 1.0_rp / (CN * CJ - CO * CI)
            imp%GA(row) = (CL * CJ - CO * CG) * mult
            imp%GB(row) = (CM * CJ - CO * CH) * mult
            imp%GC(row) = (-CP * CJ + CO * CK) * mult
            mult = 1.0_rp / (CO * CI - CN * CJ)
            imp%GD(row) = (CL * CI - CN * CG) * mult
            imp%GE(row) = (CM * CI - CN * CH) * mult
            imp%GF(row) = (-CP * CI + CN * CK) * mult
#ifdef DEBUG_L1
            if (isnan(imp%GA(row))) then
                print *, "NAN(imp%GA)", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                mult = 1.0_rp / (CN * CJ - CO * CI)
                print *, "--details:(", CL, "*", CJ, "-", CO, "*", CG, ")*", mult
                read(*,*)
            end if
            if (isnan(imp%GB(row))) then
                print *, "NAN(imp%GB)", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                mult = 1.0_rp / (CN * CJ - CO * CI)
                print *, "--details:(", CM, "*", CJ, "-", CO, "*", CH, ")*", mult
                read(*,*)
            end if
            if (isnan(imp%GC(row))) then
                print *, "NAN(imp%GC)", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                mult = 1.0_rp / (CN * CJ - CO * CI)
                print *, "--details:(", -CP, "*", CJ, "-", CO, "*", CK, ")*", mult
                read(*,*)
            end if
            if (isnan(imp%GD(row))) then
                print *, "NAN(imp%GD)", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                mult = 1.0_rp / (CO * CI - CN * CJ)
                print *, "--details:(", CL, "*", CI, "-", CN, "*", CG, ")*", mult
                read(*,*)
            end if
            if (isnan(imp%GE(row))) then
                print *, "NAN(imp%GE)", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                mult = 1.0_rp / (CO * CI - CN * CJ)
                print *, "--details:(", CM, "*", CI, "-", CN, "*", CG, ")*", mult
                read(*,*)
            end if
            if (isnan(imp%GF(row))) then
                print *, "NAN(imp%GF)", row, iseg, ie, msh%seg(iseg)%first_cs
                print *, "--cs_range:", msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
                mult = 1.0_rp / (CO * CI - CN * CJ)
                print *, "--details:(", -CP, "*", CI, "-", CN, "*", CG, ")*", mult
                read(*,*)
            end if
#endif      
#ifdef DEBUG_L3        
            print *, "preissmann_double_sweep_LPI:", iseg, ie, row
            print *, "MC:", CG, CH, CI, CJ, CK
            print *, "Q1:", CB1, CB2, CB3
            print *, "Q2:", CC1, CC4, CC5
            print *, "Q3:", CD1, CD4, CD5
            print *, "Q4:", CH1, CH4, CH5
            print *, "Q5:", CE2, CE3
            print *, "Q6:", CF1, CF2, CF3, CF4, CF5
            print *, "Q*:", CL, CM, CN, CO, CP
            print *, "GA-C:", imp%GA(row), imp%GB(row), imp%GC(row)
            print *, "GD-F:", imp%GD(row), imp%GE(row), imp%GF(row)
            !read(*,*)
#endif
            
            ! Update state variables for left cross-section
            zL = zR
            qL = qR
            aL = aR
            pL = pR
            wL = wR
            uL = uR
            debL = debR
            ddebdzL = ddebdzR
            betaL = betaR
            FrL = FrR
            ddaL = ddaR
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
            acL(:) = acR(:)
            pcL(:) = pcR(:)
            wcL(:) = wcR(:)
            kcL(:) = kcR(:)
            rhcL(:) = rhcR(:)
            dpcdzL(:) = dpcdzR(:)
            AdL = AdR
#endif
            kL = kR
            rhL = rhR
            dpdzL = dpdzR
            dkdzL = dkdzR
        end do
        
    end do
    
    ! Compute condensed transfer coefficients
#ifdef DEBUG_L3
    print *, "Compute condensed transfer relations"
#endif
    do iseg = 1, msh%nseg
    
#ifdef OLD_VERSION
        ! Compute condensed relations x1 = TA.xn + TB
        row = imp%seg_offsets(iseg) + msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs
        imp%TA1(iseg) = imp%GA(row)
        imp%TA2(iseg) = imp%GB(row)
        imp%TA3(iseg) = imp%GD(row)
        imp%TA4(iseg) = imp%GE(row)
        imp%TB1(iseg) = imp%GC(row)
        imp%TB2(iseg) = imp%GF(row)
        do ie = msh%seg(iseg)%last_cs - 2, msh%seg(iseg)%first_cs, -1
        
            row = row - 1
            imp%TA1(iseg) = imp%GA(row) * imp%TA1(iseg) + imp%GB(row) * imp%TA3(iseg)
            imp%TA2(iseg) = imp%GA(row) * imp%TA2(iseg) + imp%GB(row) * imp%TA4(iseg)
            imp%TA3(iseg) = imp%GD(row) * imp%TA1(iseg) + imp%GE(row) * imp%TA3(iseg)
            imp%TA4(iseg) = imp%GD(row) * imp%TA2(iseg) + imp%GE(row) * imp%TA4(iseg)
            imp%TB1(iseg) = imp%GA(row) * imp%TB1(iseg) + imp%GB(row) * imp%TB2(iseg) + imp%GC(row)
            imp%TB2(iseg) = imp%GD(row) * imp%TB1(iseg) + imp%GE(row) * imp%TB2(iseg) + imp%GF(row)

        end do
#else
        row = imp%seg_offsets(iseg) + 1
        mult = 1.0 / imp%GA(row)
        imp%TA1(iseg) = -mult
        imp%TA2(iseg) = imp%GB(row) * mult
        imp%TB1(iseg) = -imp%GC(row) * mult
        imp%TA3(iseg) = imp%GD(row) * imp%TA1(iseg)
        imp%TA4(iseg) = imp%GD(row) * imp%TA2(iseg) - imp%GE(row)
        imp%TB2(iseg) = imp%GF(row) + imp%GD(row) * imp%TB1(iseg)
        do ie = msh%seg(iseg)%first_cs+1, msh%seg(iseg)%last_cs-1
            row = row + 1
            
            ! Compute multiplier mult = 1 / (GA(k) + TA2(k-1) * GD(k))
            mult = 1.0_rp / (imp%GA(row) + imp%TA2(iseg) * imp%GD(row))
            
            ! Update TB1(k) = (TB1(k-1) - GC(k) + TA2(k-1) * GF(k)) / (GA(k) + TA2(k-1) * GD(k))
            imp%TB1(iseg) = (imp%TB1(iseg) - (imp%GC(row) + imp%GF(row) * imp%TA2(iseg))) * mult
            
            ! Update TA1(k) = TA1(k-1) / (GA(k) + TA2(k-1) * GD(k))
            imp%TA1(iseg) = imp%TA1(iseg) * mult
            
            ! Update TA2(k) = (GB(k) + TA2(k-1) * GE(k)) / (GA(k) + TA2(k-1) * GD(k))
            imp%TA2(iseg) = (imp%GB(row) + imp%GE(row) * imp%TA2(iseg)) * mult
            
            ! Update TB2(k) = TB2(k-1) - TA4(k-1) * (GD(k) * TB1(k) + GF(k))
            imp%TB2(iseg) = imp%TB2(iseg) - imp%TA4(iseg) * (imp%GD(row) * imp%TB1(iseg) + imp%GF(row))
            
            ! Update TA3(k) = (TA3(k-1) - TA4(k-1) * TA1(k) * GD(k))
            imp%TA3(iseg) = imp%TA3(iseg) - imp%TA4(iseg) * imp%TA1(iseg) * imp%GD(row)
            
            ! Update TA4(k) = TA4(k-1) * (GE(k) - GD(k) * TA2(k))
            imp%TA4(iseg) = imp%TA4(iseg) * (imp%GE(row) - imp%GD(row) * imp%TA2(iseg))
            
#ifdef DEBUG_L1
            if (isnan(imp%TB1(iseg))) then
                print *, "NAN(TB1),", imp%TB1(iseg), iseg, ie
                print *, "--details:(", imp%GC(row),  "+", imp%GF(row), "*", imp%TA2(iseg), ")*", mult
                call write_dof_records(0, msh, dof)
                read(*,*)
            end if
            ! Update TA1(k) = TA1(k-1) / (GA(k) + TA2(k-1) * GD(k))
            if (isnan(imp%TA1(iseg))) then
                print *, "NAN(TA1),", imp%TA1(iseg), iseg, ie
                print *, "--details:", mult
                read(*,*)
            end if
            ! Update TA2(k) = (GB(k) + TA2(k-1) * GE(k)) / (GA(k) + TA2(k-1) * GD(k))
            if (isnan(imp%TA2(iseg))) then
                print *, "NAN(TA2),", imp%TA2(iseg), iseg, ie
                print *, "--details:(", imp%GB(row), "+", imp%GE(row), "*", imp%TA2(iseg), ")*", mult
                read(*,*)
            end if
            
            if (isnan(imp%TB2(iseg))) then
                print *, "NAN(TB2),", imp%TB2(iseg), iseg, ie
                print *, "--details:", imp%TA4(iseg), "*", imp%GD(row), "*", imp%TB1(iseg), "+", imp%GF(row)
                read(*,*)
            end if
            if (isnan(imp%TA3(iseg))) then
                print *, "NAN(TA3),", imp%TA3(iseg), iseg, ie
                print *, "--details:", imp%TA4(iseg), "*", imp%TA1(iseg), "*", imp%GD(row)
                read(*,*)
            end if
            if (isnan(imp%TA4(iseg))) then
                print *, "NAN(TA4),", imp%TA4(iseg), iseg, ie
                print *, "--details:", imp%GE(row), "-", imp%GD(row), "*", imp%TA2(iseg)
                read(*,*)
            end if
#endif            
        end do
#endif
        
    end do
    
    ! Compute reduced system
    off = 1
    row = 1
    do iseg = 1, msh%nseg
    
        ! Upstream boundary condition
        if (msh%seg(iseg)%us_bc > 0) then
        
            ibc = msh%seg(iseg)%us_bc
            
            if (mdl%bc(ibc)%id == "discharge") then
            
                imp%ANZ(off) = 1.0_rp
                imp%RHS(row) = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, mdl%tc) - dof%Q(msh%seg(iseg)%first_cs)
                off = off + 1
                
            end if
            
        else
        
            imp%RHS(row) = dof%q(msh%seg(iseg)%first_cs)
            do ius_seg = 1, size(msh%seg(iseg)%us_seg)
                us_seg = msh%seg(iseg)%us_seg(ius_seg)
                imp%ANZ(off) = 1.0_rp
                off = off + 1
                imp%RHS(row) = imp%RHS(row) - dof%q(msh%seg(us_seg)%last_cs)
            end do
            imp%ANZ(off) = -1.0_rp
            off = off + 1
            
        end if
        row = row + 1
        
        ! Condensed transfer relations
#ifdef OLD_VERSION
        imp%ANZ(off) = 1.0_rp
        imp%ANZ(off+1) = 0.0_rp
        imp%ANZ(off+2) = -imp%TA1(iseg)
        imp%ANZ(off+3) = -imp%TA2(iseg)
        imp%RHS(row) = imp%TB1(iseg)
        off = off + 4
        row = row + 1
        imp%ANZ(off) = 0.0_rp
        imp%ANZ(off+1) = 1.0_rp
        imp%ANZ(off+2) = -imp%TA3(iseg)
        imp%ANZ(off+3) = -imp%TA4(iseg)
        imp%RHS(row) = imp%TB2(iseg)
        off = off + 4
        row = row + 1
#else
        imp%ANZ(off) = imp%TA1(iseg)
        imp%ANZ(off+1) = 0.0_rp
        imp%ANZ(off+2) = 1.0_rp
        imp%ANZ(off+3) = imp%TA2(iseg)
        imp%RHS(row) = imp%TB1(iseg)
        off = off + 4
        row = row + 1
        imp%ANZ(off) = imp%TA3(iseg)
        imp%ANZ(off+1) = 1.0_rp
        imp%ANZ(off+2) = 0.0_rp
        imp%ANZ(off+3) = imp%TA4(iseg)
        imp%RHS(row) = imp%TB2(iseg)
#ifdef DEBUG_L1
        if (isNan(imp%ANZ(off))) then
            print *, "imp%ANZ(off)=,", imp%ANZ(off), iseg
            print *, "--TA3,", imp%TA3(iseg)
            read(*,*)
        end if
        if (isNan(imp%ANZ(off+1))) then
            print *, "imp%ANZ(off+1)=,", imp%ANZ(off+3), iseg
!             print *, "--TA4,", imp%TA4(iseg)
            read(*,*)
        end if
        if (isNan(imp%ANZ(off+2))) then
            print *, "imp%ANZ(off+2)=,", imp%ANZ(off+3), iseg
!             print *, "--TA4,", imp%TA4(iseg)
            read(*,*)
        end if
        if (isNan(imp%ANZ(off+3))) then
            print *, "imp%ANZ(off+3)=,", imp%ANZ(off+3), iseg
            print *, "--TA4,", imp%TA4(iseg)
            read(*,*)
        end if
#endif
        off = off + 4
        row = row + 1
#endif
        
        ! Downstream boundary condition
        if (msh%seg(iseg)%ds_bc > 0) then
            ibc = msh%seg(iseg)%ds_bc
            if (mdl%bc(ibc)%id == "elevation") then
                imp%ANZ(off) = 1.0_rp
                zL = dof%h(msh%seg(iseg)%last_cs) + msh%cs(msh%seg(iseg)%last_cs)%bathy
                imp%RHS(row) = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, mdl%tc) - zL
!                 print *, "DOWNSTREAM:", zl, imp%RHS(row), row, imp%ANZ(off)
!                 print *, "tc:", mdl%tc, msh%seg(iseg)%last_cs
!                 print *, "ts%t:", mdl%bc(ibc)%ts%t
!                 print *, "ts%y:", mdl%bc(ibc)%ts%y
!                 read(*,*)
                off = off + 1
            else if (mdl%bc(ibc)%id == "normal_depth") then
                ics = msh%seg(iseg)%last_cs
                So = msh%cs(ics)%slope
                if (msh%strickler_type_code == strickler_type_average_conveyance) then
#ifdef AVERAGE_CONVEYANCE
                    call areas_compound_channel(msh%cs(ie), dof%h(ie), acL)
                    call perimeters_compound_channel(msh%cs(ie), dof%h(ie), pcL)
                    call widths_compound_channel(msh%cs(ie), dof%h(ie), wcL)
                    kcL(:) = msh%cs(ie)%strickler_params(:)
                    call dpdz_compound_channel(msh%cs(ie), dof%h(ie), dpcdzL)
                    do i = 1, 3
                        dpcdzL(i) = min(100.0_rp, dpcdzL(i))
                        if (pcL(i) > 1.e-12_rp) then
                            rhcL(i) = acL(i) / pcL(i)
                        else
                            rhcL(i) = 0.d0
                        end if
                    end do
                    ddebdzL = d1p3 * (kcL(1) * rhcL(1)**d2p3 * (5.0_rp * wcL(1) - 2.0_rp * rhcL(1) * dpcdzL(1)) + &
                                      kcL(2) * rhcL(2)**d2p3 * (5.0_rp * wcL(2) - 2.0_rp * rhcL(2) * dpcdzL(2)) + &
                                      kcL(3) * rhcL(3)**d2p3 * (5.0_rp * wcL(3) - 2.0_rp * rhcL(3) * dpcdzL(3)))
#endif
#ifdef DEBORD_FORMULA
                else if (msh%strickler_type_code == strickler_type_debord) then

                    ! Retrieve compound state variables
                    kcL(:) = msh%cs(ie)%strickler_params(:)
                    call areas_compound_channel(msh%cs(ie), dof%h(ie), acL)
                    call perimeters_compound_channel(msh%cs(ie), dof%h(ie), pcL)
                    call widths_compound_channel(msh%cs(ie), dof%h(ie), wcL)
                    call dpdz_compound_channel(msh%cs(ie), dof%h(ie), dpcdzL)
                    
                    ! Compute variables for floodplain (composition of right and left overbanks)
                    if (pcL(1) + pcL(3) > zerom) then
                        kcL(1) = (kcL(1) * pcL(1) + kcL(3) * pcL(3)) / (pcL(1) + pcL(3))
                    else
                        kcL(1) = 0.5_rp * (kcL(1) + kcL(3))
                    endif
                    pcL(1) = pcL(1) + pcL(3)
                    acL(1) = acL(1) + acL(3)
                    wcL(1) = wcL(1) + wcL(3)
                    dpcdzL(1) = dpcdzL(1) + dpcdzL(3)

                    ! Compute hydraulic radiuses and limit dPdZ
                    do i = 1, 2
                        dpcdzL(i) = min(mdl%dpdz_cutoff, dpcdzL(i))
                        if (pcL(i) > tinym) then
                            rhcL(i) = acL(i) / pcL(i)
                        else
                            rhcL(i) = 0.d0
                        end if
                    end do
                    
                    ! Compute A0 and A Debord
                    A0d = 0.9_rp * (kcL(1) / kcL(2))**d1p6
                    if (rhcL(1) > tinym) then
                        r = rhcL(1) / rhcL(2)
                        if (r > 0.3_rp) then
                            Ad = A0d
                        else
                            Ad = demi * (1.0_rp - A0d) * cos(pi * r / 0.3_rp) + demi * (1.0_rp + A0d)
                        end if
                    else
                        Ad = 1.0_rp
                    end if
                    AdL = Ad
                    
                    ! Update Stricklers before conveyance computation
                    if (acL(2) > tinym) then
                        kcL(1) = sqrt(1.0_rp + acL(1) / acL(2) * (1.0_rp - Ad**2)) * kcL(1)
                        kcL(2) = Ad * kcL(2)
                    end if
                    
                    ! Compute total gradient of conveyance w.r.t. to height
                    ddebdzL = d1p3 * (kcL(1) * rhcL(1)**d2p3 * (5.0_rp * wcL(1) - 2.0_rp * rhcL(1) * dpcdzL(1)) + &
                                      kcL(2) * rhcL(2)**d2p3 * (5.0_rp * wcL(2) - 2.0_rp * rhcL(2) * dpcdzL(2)))
#endif
                else
                    call area(msh%cs(ics), dof%h(ics), aL)
                    call perimeter(msh%cs(ics), dof%h(ics), pL)
                    call dpdz(msh%cs(ics), dof%h(ics), dpdzL)
                    dpdzL = min(100.0_rp, dpdzL)
                    call dkdz(msh%cs(ie), dof%h(ie), dkdzL)
                    call strickler(msh%cs(ics), dof%h(ics), kL)
                    rhL = aL / pL
                    ddebdzL = dkdzL * aL * rhL**d2p3 + &
                            d1p3 * (kL * rhL**d2p3) * (5.0_rp * wL - 2.0_rp * rhL * dpdzL)
                end if
                imp%ANZ(off) = 1.0_rp
                imp%ANZ(off+1) = -ddebdzL * sqrt(So)
#ifdef DEBUG_L1
                if (isNan(imp%ANZ(off+1))) then
                    print *, "BC[normal_depth]:imp%ANZ(off+1)=,", imp%ANZ(off+1)
                    print *, "--kL,", kL
                    print *, "--rhL,", rhL
                    print *, "--wL,", wL
                    print *, "--dpdzL,", dpdzL
                    read(*,*)
                end if
#endif
                imp%RHS(row) = 0.0_rp
                off = off + 2
            else if (mdl%bc(ibc)%id == "rating_curve") then
                ics = msh%seg(iseg)%last_cs
                zL = dof%h(ics) + msh%cs(ics)%bathy- 0.1_rp
                zR = dof%h(ics) + msh%cs(ics)%bathy+ 0.1_rp
                qL = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, zL)
                qR = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, zR)
                imp%ANZ(off) = 1.0_rp
                imp%ANZ(off+1) = -(qR - qL) / (zR - zL)
                imp%RHS(row) = 0.0_rp
                off = off + 2
            end if
        else
            imp%ANZ(off) = 1.0_rp
            imp%ANZ(off+1) = -1.0_rp
!             imp%RHS(row) = 0.0_rp
            ds_seg = msh%seg(iseg)%ds_seg
            imp%RHS(row) = dof%h(msh%seg(ds_seg)%first_cs) + msh%cs(msh%seg(ds_seg)%first_cs)%bathy - &
                           (dof%h(msh%seg(iseg)%last_cs) + msh%cs(msh%seg(iseg)%last_cs)%bathy)
            off = off + 2
        end if
        row = row + 1

    end do

    ! Solve reduced system
    call solve_system_mumps(0, imp)
    
    ! Compute solution using double sweep on each segment
    do iseg = 1, msh%nseg
    
        ics = msh%seg(iseg)%first_cs
        row = imp%seg_offsets(iseg) + 1
        imp%CR(row) = 1.0_rp
        imp%CS(row) = 0.0_rp
        imp%CT(row) = imp%RHS(4 * (iseg - 1) + 1)
        do ics = msh%seg(iseg)%first_cs+1, msh%seg(iseg)%last_cs
            imp%CR(row+1) = imp%CR(row) * imp%GA(row) + imp%CS(row) * imp%GD(row)
            imp%CS(row+1) = imp%CR(row) * imp%GB(row) + imp%CS(row) * imp%GE(row)
            imp%CT(row+1) = imp%CT(row) - (imp%CR(row) * imp%GC(row) + imp%CS(row) * imp%GF(row))
            if (abs(imp%CR(row+1)) > 1e-02) then
              imp%CS(row+1) = imp%CS(row+1) / imp%CR(row+1)
              imp%CT(row+1) = imp%CT(row+1) / imp%CR(row+1)
              imp%CR(row+1) = 1.0_rp
            end if
            row = row + 1
        end do
        
        ics = msh%seg(iseg)%last_cs
        row = imp%seg_offsets(iseg) + msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs
        DQ(ics) = imp%RHS(4 * (iseg - 1) + 3)
        DZ(ics) = imp%RHS(4 * (iseg - 1) + 4)
        dof%h(ics) = dof%h(ics) + DZ(ics)
        if (dof%h(ics) < mdl%heps) then
            dof%h(ics) = mdl%heps
#ifndef CPP_ADJ
            mdl%warning_counters(1) = mdl%warning_counters(1) + 1
#endif
        end if
! #ifndef CPP_ADJ
!         if (msh%cs(ics)%ob_levels(1) > 0 .or. msh%cs(ics)%ob_levels(2) > 0) then
!             if (dof%h(ics) > msh%cs(ics)%level_heights(minval(msh%cs(ics)%ob_levels)) - msh%cs(ics)%bathy) then
!                 print *, "overbank flow:", ics
!             end if
!         end if
! #endif
        dof%A(ics) = htoA(msh%cs(ics), dof%h(ics))
        dof%Q(ics) = max(mdl%qeps, dof%Q(ics) + DQ(ics))
        if (isnan(dof%h(ics)) .or. isnan(dof%Q(ics))) then
#ifndef CPP_ADJ
!             print *, "[ SIMULATION FAILED ] NaN depth or discharge"
#endif
            status = 21
            open(199, file="error.txt")
            write(199, '(2(A,I6))') "NaN dof detected on segment ", iseg-1, ", cross-section ", ics
            close(199)

#ifdef DEBUG_L1
            print *, "NAN,", ics, ", h=", dof%h(ics), ", q=", dof%q(ics), DZ(ics), DQ(ics)
            read(*,*)
            print *, imp%ANZ
            read(*,*)
#endif
            return
        end if
        ics = msh%seg(iseg)%first_cs
        do ics = msh%seg(iseg)%last_cs-1, msh%seg(iseg)%first_cs, -1
            DZ(ics) = imp%GD(row) * DQ(ics+1) + imp%GE(row) * DZ(ics+1) + imp%GF(row)
            DQ(ics) = (imp%CT(row) - imp%CS(row) * DZ(ics)) / imp%CR(row)
            dof%h(ics) = dof%h(ics) + DZ(ics)
            if (dof%h(ics) < mdl%heps) then
                dof%h(ics) = mdl%heps
#ifndef CPP_ADJ
                mdl%warning_counters(1) = mdl%warning_counters(1) + 1
#endif
            end if
! #ifndef CPP_ADJ
!             if (msh%cs(ics)%ob_levels(1) > 0 .or. msh%cs(ics)%ob_levels(2) > 0) then
!                 if (msh%cs(ics)%ob_levels(1) > 0 .and. &
!                     dof%h(ics) > msh%cs(ics)%level_heights(msh%cs(ics)%ob_levels(1)) - msh%cs(ics)%bathy) then
!                     print *, "overbank flow:", ics, dof%h(ics)+msh%cs(ics)%bathy
!                 else if (msh%cs(ics)%ob_levels(2) > 0 .and. &
!                     dof%h(ics) > msh%cs(ics)%level_heights(msh%cs(ics)%ob_levels(2)) - msh%cs(ics)%bathy) then
!                     print *, "overbank flow:", ics, dof%h(ics)+msh%cs(ics)%bathy
!                 end if
!             end if
! #endif
            dof%A(ics) = htoA(msh%cs(ics), dof%h(ics))
            dof%Q(ics) = max(mdl%qeps, dof%Q(ics) + DQ(ics))
            if (isnan(dof%h(ics)) .or. isnan(dof%Q(ics))) then
! #ifndef CPP_ADJ
!                 print *, "[ SIMULATION FAILED ] NaN depth or discharge"
! #endif
                status = 21
                open(199, file="error.txt")
                write(199, '(2(A,I6))') "NaN dof detected on segment ", iseg-1, ", cross-section ", ics
                close(199)

#ifdef DEBUG_L1
                print *, "NAN,", ics, ", h=", dof%h(ics), ", q=", dof%q(ics), DZ(ics), DQ(ics)
                read(*,*)
#endif
                return
            end if
            row = row - 1
        
        end do
        
    end do
#ifdef DEBUG_L1
    has_too_big_depth = .false.
    do iseg = 1, msh%nseg
        do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
            if (isnan(dof%h(ics))) then
                print *, "FINAL CHECK;h is NaN:", iseg, ics
                read(*,*)
            end if
            if (dof%h(ics) - (mdl%msh%cs(ics)%level_heights(1) - mdl%msh%cs(ics)%bathy)  > 100.0) then
                has_too_big_depth = .true.
                print *, "FINAL CHECK;h is TOO BIG:", iseg, ics, dof%h(ics)
    !             read(*,*)
            end if
        end do
    end do
    do iseg = 1, size(msh%seg)
        do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
            if (dof%h(ics) < 0.01) then
                print *, "FINAL CHECK;h is TOO SMALL:", ics, dof%h(ics)
    !             read(*,*)
            end if
        end do
    end do
    if (has_too_big_depth) then
        print *, "write debug_preissmann.dat at t=", mdl%tc
        open(999, file="debug_preissmann.dat")
        do ics = 3, size(mdl%msh%cs)-2
            write(999, *) mdl%msh%cs(ics)%x, mdl%msh%cs(ics)%bathy, mdl%msh%cs(ics)%bathy + mdl%dof%h(ics), mdl%msh%cs(ics)%bathy + hprec(ics), dof%q(ics), qprec(ics)
        end do
        close(999)
        read(*,*)
    end if
#endif
#ifndef CPP_ADJ
!     do ics = 1, size(msh%cs)
!         if (.not. IEEE_IS_FINITE(dof%h(ics)) .or. &
!             .not. IEEE_IS_FINITE(dof%Q(ics)) .or. &
!             dof%h(ics) - (mdl%msh%cs(ics)%level_heights(1) - mdl%msh%cs(ics)%bathy)  > 100.0) then
!             status = -1
!         end if
!     end do
#endif

!     print *, "PAUSE..."
!     read(*,*)
   
end subroutine preissmann_timestep





subroutine write_dof_records(ival, msh, dof)
    use m_common
    use m_mesh
    use m_sw_mono, only: Unknowns
    use m_numeric
#ifdef DEBUG_L1
    use, intrinsic :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param msh Mesh of the model
    integer(ip), intent(in) :: ival
    !> @param msh Mesh of the model
    TYPE(Mesh), intent(inout) :: msh
    !> @param dof Unknowns of the model
    TYPE(Unknowns), intent(inout) :: dof
    !=LOCAL VARIABLES==================================================================================================!
    integer(ip) :: i
    
    open(10, file="dof_record.dat")
    do i = 3, msh%ncs-2
        write(10, '(E12.6,10(1X,E12.6))') msh%cs(i)%x, msh%cs(i)%bathy, msh%cs(i)%level_heights(1), dof%h(i), dof%a(i), dof%q(i)
    end do
    close(10)
    
end subroutine
