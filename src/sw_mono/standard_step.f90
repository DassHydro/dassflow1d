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
!> @file standard_step.f90
!> @brief This file includes the routine that solves the Steady State Shallow Water Equation using standard step method

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compute backwater curve
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

#ifndef CPP_ADJ
! #define DEBUG_L1
! #define DEBUG_L2
! #define DEBUG_SAVE_PROFILE
#endif 

!> @brief Solve the Steady State Shallow Water Equation using standard step method
subroutine standard_step(mdl, eps, itermax)
    use m_mesh
    use m_sw_mono, only: Model
    use m_numeric
!     use m_user_data
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Model), intent(inout) :: mdl
    real(rp), intent(in) :: eps
    integer(ip), intent(in) :: itermax
    !=LOCAL VARIABLES==================================================================================================!
    logical :: is_upstream_seg
    ! Index of downstream segment
    integer(ip) :: ds_seg
    ! Index of first cross-section
    integer(ip) :: first_cs
    ! Iterator
    integer(ip) :: i
    ! Index of boundary condition
    integer(ip) :: ibc
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of inflow condition
    integer(ip) :: iic
    ! Index of segment
    integer(ip) :: iseg
    ! Iterations counter
    integer(ip) :: iter
    ! Index of last cross-section
    integer(ip) :: last_cs
    ! Index of segment
    integer(ip) :: seg_stack_size
    ! Index of upstream segment
    integer(ip) :: us_seg
    ! Stack of segments
    integer(ip), dimension(mdl%msh%nseg) :: seg_stack
    ! State variables for boundary cross-section
    real(rp) :: aB, kB, pB, rhB, wB, dpdzB, dkdzB
    ! Increment in depth
    real(rp) :: dh
    ! Spacing
    real(rp) :: dx
    ! Functional and gradient for Newton loop
    real(rp) :: F, dF
    ! Multiplier
    real(rp) :: mult
    ! Discharge
    real(rp) :: q
    ! Lateral inflow
    real(rp) :: qlat
    ! Heights for secant search
    real(rp) :: z1, z2, zc
    ! Bathymetry slope
    real(rp) :: So
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
    ! State variables for boundary cross-section for compound channels
    real(rp), dimension(3) :: acB, kcB, pcB, rhsB, wcB, dpcdzB
#endif
#ifdef DEBUG_L1
    logical :: has_too_big_depth
    real(rp) :: max_h0
#endif

    ! Set interior discharge values to zero for all segments
    do iseg = 1, mdl%msh%nseg
    
        mdl%dof%q(mdl%msh%seg(iseg)%first_cs:mdl%msh%seg(iseg)%last_cs) = 0.0_rp
        
    end do
    
    ! Compute discharges
    do iseg = 1, mdl%msh%nseg
    
        is_upstream_seg = .false.
    
        if (.not.allocated(mdl%msh%seg(iseg)%us_seg)) then
        
            is_upstream_seg = .true.
            
        else if (mdl%msh%seg(iseg)%us_seg(1) < 1) then
        
            is_upstream_seg = .true.

        end if
        
        if (is_upstream_seg) then
        
            ibc = mdl%msh%seg(iseg)%us_bc
            
            if (mdl%bc(ibc)%id == "discharge") then
            
                q = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, mdl%tc)
                ! print *, "standard step: q=", q
                
            end if
        
            mdl%dof%q(mdl%msh%seg(iseg)%first_cs-2:mdl%msh%seg(iseg)%last_cs+2) = q
            
            ! Append discharge for downstream segments
            ds_seg = mdl%msh%seg(iseg)%ds_seg
            
            do while(ds_seg > 0)
            
                q = mdl%dof%q(mdl%msh%seg(ds_seg)%first_cs) + mdl%dof%q(mdl%msh%seg(iseg)%first_cs-1)
                mdl%dof%q(mdl%msh%seg(ds_seg)%first_cs:mdl%msh%seg(ds_seg)%last_cs) = q
                ds_seg = mdl%msh%seg(ds_seg)%ds_seg
            
            end do
            
        end if
        
    end do
    
    ! Apply inflow discharges
    if (allocated(mdl%ic) .and. mdl%nic>0) then
    
        do iic = 1, size(mdl%ic)
        
            ! Update discharges for segment where inflow applies
            iseg = mdl%ic(iic)%iseg
            ics = mdl%ic(iic)%ie+1
            mult =  1.0 / mdl%ic(iic)%ledge
            qlat = linear_interp(mdl%ic(iic)%ts%t, mdl%ic(iic)%ts%y, mdl%tc) * mdl%msh%cs(ics)%deltademi
            q = mdl%dof%q(ics) + qlat * mult
            mdl%dof%q(ics:mdl%msh%seg(iseg)%last_cs) = q
                    
            ! Update discharges for downstream segments
            ds_seg = mdl%msh%seg(iseg)%ds_seg
            
            do while(ds_seg > 0)
            
                do ics = mdl%msh%seg(ds_seg)%first_cs, mdl%msh%seg(ds_seg)%last_cs
                    mdl%dof%q(ics) = mdl%dof%q(ics) + qlat * mult
                end do
                ds_seg = mdl%msh%seg(ds_seg)%ds_seg
            
            end do
            
        end do
    end if

    ! Retrieve downstream segment
    ds_seg = 0
    do iseg = 1, mdl%msh%nseg
    
        if (mdl%msh%seg(iseg)%ds_seg <= 0) then
        
#ifndef CPP_ADJ
            ! Check that a previous downstream segment has not been found
            if (ds_seg /= 0) then
                call abort_solver("Multiple downstream segments in mesh")
            end if
#endif
            
            ! Set downstream segment index
            ds_seg = iseg
            
        end if
    
    end do
    
    ! Compute boundary elevation
    ibc = mdl%msh%seg(ds_seg)%ds_bc
    ics = mdl%msh%seg(ds_seg)%last_cs
    if (mdl%bc(ibc)%id == "elevation") then
    
        if (allocated(mdl%bc(ibc)%ts%t)) then
            mdl%dof%h(ics) = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, mdl%tc) - mdl%msh%cs(ics)%bathy
            call update_level(mdl%msh%cs(ics), mdl%dof%h(ics))
        else
            mdl%dof%h(ics) = mdl%bc(ibc)%ts%y(1) - mdl%msh%cs(ics)%bathy
            call update_level(mdl%msh%cs(ics), mdl%dof%h(ics))
        end if
        mdl%dof%A(ics) = HtoA(mdl%msh%cs(ics), mdl%dof%h(ics))
        
#ifndef CPP_ADJ
        ! Check for NaN value(s)
        if (isnan(mdl%dof%h(ics)) .or. isnan(mdl%dof%A(ics))) then
            mdl%status = 10
            return
        end if
#endif
#ifdef DEBUG_L1
        print *, "DN", ics, mdl%dof%h(ics), mdl%dof%a(ics), mdl%dof%q(ics)
#endif
        
    else if (mdl%bc(ibc)%id == "normal_depth") then
    
!         print *, "NORMAL_DEPTH:BEGIN"
        mdl%dof%h(ics) = 1.0
        mdl%dof%A(ics) = htoA_noupdate(mdl%msh%cs(ics), mdl%dof%h(ics))
        dh = 1.0
        dx = mdl%msh%cs(ics-1)%deltademi
!         So = max(1e-8, (mdl%msh%cs(ics-1)%bathy - mdl%msh%cs(ics)%bathy) / dx)
        So = max(1e-8, mdl%msh%cs(ics)%slope)
        iter = 0
#ifdef DEBUG_L1
        print *, "NORMAL_DEPTH:So=", So, ics, ics-1
#endif
!         print *, "NORMAL_DEPTH:First pass..."
        do while (abs(dh) > 1e-3 .and. iter < 100)
        
            call update_level(mdl%msh%cs(ics), mdl%dof%h(ics))
            
            if (mdl%msh%strickler_type_code == strickler_type_average_conveyance .or. &
                mdl%msh%strickler_type_code == strickler_type_debord) then
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
                call perimeters_compound_channel(mdl%msh%cs(ics), mdl%dof%h(ics), pcB)
                call areas_compound_channel(mdl%msh%cs(ics), mdl%dof%h(ics), acB)
                call widths_compound_channel(mdl%msh%cs(ics), mdl%dof%h(ics), wcB)
                call dPdZ_compound_channel(mdl%msh%cs(ics), mdl%dof%h(ics), dpcdzB)
                kcB(:) = mdl%msh%cs(ics)%strickler_params(:)
                rhsB(:) = acB(:) / max(1.e-12_rp, pcB(:))
                dF = (rhsB(1)**d2p3 * kcB(1) * (d5p3 * wcB(1) - d2p3 * rhsB(1) * dpcdzB(1)) + &
                      rhsB(2)**d2p3 * kcB(2) * (d5p3 * wcB(2) - d2p3 * rhsB(2) * dpcdzB(2)) + &
                      rhsB(3)**d2p3 * kcB(3) * (d5p3 * wcB(3) - d2p3 * rhsB(3) * dpcdzB(3))) * sqrt(So)
                F = (kcB(1) * acB(1) * rhsB(1)**d2p3 + &
                     kcB(2) * acB(2) * rhsB(2)**d2p3 + &
                     kcB(3) * acB(3) * rhsB(3)**d2p3) * sqrt(So) - mdl%dof%q(ics)
#ifdef DEBUG_L1
                print *, "--Hn", mdl%dof%h(ics)
                print *, "--Q=", mdl%dof%q(ics)
                print *, "--A=", acB
                print *, "--P=", pcB
                print *, "--W=", wcB
                print *, "--Rh=", rhsB
                print *, "--So=", So
                print *, "--Qn=", F
                print *, "--dh=", dh, abs(dh), abs(dh) > 1e-3
                print *, "--F=", F
                print *, "--dF=", dF, dKdZB, dPdZB
#endif
#endif
            else
                call strickler(mdl%msh%cs(ics), mdl%dof%h(ics), kB)
                call perimeter(mdl%msh%cs(ics), mdl%dof%h(ics), pB)
                call area(mdl%msh%cs(ics), mdl%dof%h(ics), aB)
                call width(mdl%msh%cs(ics), mdl%dof%h(ics), wB)
                call dKdZ(mdl%msh%cs(ics), mdl%dof%h(ics), dkdzB)
                call dPdZ(mdl%msh%cs(ics), mdl%dof%h(ics), dpdzB)
                dpdzB = min(100.0_rp, dpdzB)
                rhB = aB / pB
                dF = rhB**d2p3 * (dkdzB * aB + kB * (d5p3 * wB - d2p3 * rhB * dpdzB)) * sqrt(So)
                F = kB * aB * rhB**d2p3 * sqrt(So) - mdl%dof%q(ics)
#ifdef DEBUG_L1
                print *, "--Hn", mdl%dof%h(ics)
                print *, "--Q=", mdl%dof%q(ics)
                print *, "--A=", aB
                print *, "--P=", pB
                print *, "--W=", wB
                print *, "--Rh=", rhB
                print *, "--K=", kB
                print *, "--So=", So
                print *, "--D=", kB * aB * rhB**d2p3
                print *, "--Qn=", kB * aB * rhB**d2p3*sqrt(So)
#endif
            end if
            dh = - F / dF
#ifdef DEBUG_L1
            print *, "--dh=", dh, abs(dh), abs(dh) > 1e-3
            print *, "--F=", F
            print *, "--dF=", dF, dKdZB, dPdZB
            read(*,*)
#endif
            mdl%dof%h(ics) = mdl%dof%h(ics) - F / dF
            if (mdl%dof%h(ics) < mdl%heps) then
                mdl%dof%h(ics) = mdl%heps
#ifndef CPP_ADJ
                mdl%warning_counters(1) = mdl%warning_counters(1) + 1
#endif
            end if
            mdl%dof%A(ics) = HtoA(mdl%msh%cs(ics), mdl%dof%h(ics))
            
#ifndef CPP_ADJ
            ! Check for NaN value(s)
            if (isnan(mdl%dof%h(ics)) .or. isnan(mdl%dof%A(ics))) then 
                print *, "h=", mdl%dof%h(ics)
                print *, "h=", mdl%dof%A(ics)
                mdl%status = 10
                return
            end if
#endif
            iter = iter + 1
        end do
        call update_level(mdl%msh%cs(ics), mdl%dof%h(ics))
#ifdef DEBUG_L1
        print *, "DN", ics, mdl%dof%h(ics), mdl%dof%a(ics)
        read(*,*)
#endif
        
    else if (mdl%bc(ibc)%id == "rating_curve") then
    
        z1 = mdl%bc(ibc)%ts%t(1)
        z2 = mdl%bc(ibc)%ts%t(size(mdl%bc(ibc)%ts%t))
        do while(z2 - z1 > 1.e-3)
            zc = 0.5 * (z1 + z2)
            q = linear_interp(mdl%bc(ibc)%ts%t, mdl%bc(ibc)%ts%y, zc)
            if (q < mdl%dof%q(ics)) then
                z1 = zc
            else
                z2 = zc
            end if
!             print *, z1, z2, zc, q
        end do
        mdl%dof%h(ics) = 0.5 * (z1 + z2) - mdl%msh%cs(ics)%bathy
        mdl%dof%A(ics) = htoA(mdl%msh%cs(ics), mdl%dof%h(ics))
            
#ifndef CPP_ADJ
        ! Check for NaN value(s)
        if (isnan(mdl%dof%h(ics)) .or. isnan(mdl%dof%A(ics))) then
            mdl%status = 10
            return
        end if
#endif
        call update_level(mdl%msh%cs(ics), mdl%dof%h(ics))
        
    end if
    
    ! Initialise stack
    seg_stack(:) = 0
    seg_stack(1) = ds_seg
    seg_stack_size = 1
    
    ! Main loop
    do while(seg_stack_size > 0)
    
        ! Retrieve first segment in stack
        iseg = seg_stack(1)
!         print *, "iseg=", iseg
        
        ! Compute standard step on segment
        call standard_step_segment(eps, itermax, mdl, iseg)
!         if (isnan(mdl%dof%h(mdl%msh%seg(iseg)%first_cs))) then
!             print *, "FINAL CHECK(STD_STEP);h is NaN:", iseg, mdl%msh%seg(iseg)%first_cs
!             read(*,*)
!         end if
        
!         print *, "BEFORE_POP=(", seg_stack, ")", seg_stack_size
        ! Pop current segment from stack
        if (seg_stack_size > 1) then
          seg_stack(1:seg_stack_size-1) = seg_stack(2:seg_stack_size)
          seg_stack(seg_stack_size) = 0
        end if
        seg_stack_size = seg_stack_size - 1
!         print *, "AFTER_POP=(", seg_stack, ")", seg_stack_size
!         read(*,*)

        ! Set downstream elevation of upstream segments (equality of elevations at junctions)
        if (allocated(mdl%msh%seg(iseg)%us_seg)) then
          
            ! Retrieve index of first cross-section of current segment
            first_cs = mdl%msh%seg(iseg)%first_cs
            
            do i = 1, size(mdl%msh%seg(iseg)%us_seg)
                if (mdl%msh%seg(iseg)%us_seg(i) > 0) then
                
                    ! Retrieve index of current upstream segment
                    us_seg = mdl%msh%seg(iseg)%us_seg(i)
                    
                    ! Retrieve index of last cross-section of current upstream segment
                    last_cs = mdl%msh%seg(us_seg)%last_cs
                    
                    ! Compute depth of last cross-section of current upstream segment using Z_upstream = Z_downstream
                    mdl%dof%h(last_cs) = mdl%dof%h(first_cs) + mdl%msh%cs(first_cs)%bathy - mdl%msh%cs(last_cs)%bathy
                    if (mdl%dof%h(last_cs) < mdl%heps) mdl%dof%h(last_cs) = mdl%heps
#ifdef DEBUG_L1
                    print *, "UPSTREAM:", iseg, us_seg
                    print *, "- segments:", iseg, us_seg
                    print *, "- sections:", first_cs, last_cs
                    print *, "- h       :", mdl%dof%h(first_cs), mdl%dof%h(last_cs)
                    print *, "- bathy   :", mdl%msh%cs(first_cs)%bathy, mdl%msh%cs(last_cs)%bathy
                    if (mdl%dof%h(last_cs) < 0.0) then
                        read(*,*)
                    end if
#endif
                    mdl%dof%A(last_cs) = htoA(mdl%msh%cs(last_cs), mdl%dof%h(last_cs))
!                     print *, "mdl%dof%h(last_cs):", us_seg, mdl%dof%h(last_cs)
                    call update_level(mdl%msh%cs(last_cs), mdl%dof%h(last_cs))

                end if
            end do
            
        end if
!         read(*,*)
        
        ! Append stack with upstream segments
        if (allocated(mdl%msh%seg(iseg)%us_seg)) then
          
            do i = 1, size(mdl%msh%seg(iseg)%us_seg)
                if (mdl%msh%seg(iseg)%us_seg(i) > 0) then
                    seg_stack_size = seg_stack_size + 1
                    seg_stack(seg_stack_size) = mdl%msh%seg(iseg)%us_seg(i)
                end if
            end do
            
        end if
        
    end do
    
#ifdef DEBUG_L1
    max_h0 = 0.0
    do ics = 1, size(mdl%msh%cs)
        max_h0 = max(max_h0, mdl%msh%cs(ics)%level_heights(1)-mdl%msh%cs(ics)%bathy)
    end do
    print *, "MAX_H0=", max_h0

    has_too_big_depth = .false.
    do ics = 1, size(mdl%msh%cs)
        if (isnan(mdl%dof%h(ics))) then
            print *, "FINAL CHECK(STD_STEP);h is NaN:", ics
            read(*,*)
        end if
        if (mdl%dof%h(ics) > 100.0) then
            print *, "FINAL CHECK(STD_STEP);h is TOO BIG:", ics, mdl%dof%h(ics)
            print *, "--H0=", mdl%msh%cs(ics)%level_heights(1)-mdl%msh%cs(ics)%bathy
            has_too_big_depth = .true.
!             read(*,*)
        end if
    end do
    if (has_too_big_depth) then
        open(999, file="debug_standard_step.dat")
        do ics = 1, size(mdl%msh%cs)
            write(999, *) mdl%msh%cs(ics)%x, mdl%msh%cs(ics)%bathy, mdl%msh%cs(ics)%bathy + mdl%dof%h(ics)
        end do
        close(999)
        read(*,*)
    end if
#endif
#ifdef DEBUG_SAVE_PROFILE
        open(999, file="standard_step_profile.dat")
        do ics = 1, size(mdl%msh%cs)
            write(999, '(5(1X,E12.6))') mdl%msh%cs(ics)%x, mdl%msh%cs(ics)%bathy, mdl%msh%cs(ics)%bathy + mdl%dof%h(ics), &
                                        mdl%msh%cs(ics)%level_widths(1), mdl%dof%q(ics)
        end do
        close(999)
        print *, "standard_step_profile.dat UP-TO-DATE"
#endif
    
end subroutine


!> @brief Solve the Steady State Shallow Water Equation on a single segment using standard step method
subroutine standard_step_segment(eps, itermax, mdl, iseg)
    use m_mesh
    use m_sw_mono
!     use m_numeric
!     use m_user_data
    implicit none
    !=ARGUMENTS========================================================================================================!
    real(rp), intent(in) :: eps
    integer(ip), intent(in) :: itermax
    type(Model), intent(inout) :: mdl
    integer(ip), intent(in) :: iseg
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of downstream segment
    integer(ip) :: ds_seg
    ! Iterator
    integer(ip) :: i
    ! Index of boundary condition
    integer(ip) :: ibc
    ! Index of cross-section
    integer(ip) :: ics
    ! Number of iterations of convergence loop
    integer(ip) :: iter
    ! Spacing
    real(rp) :: dx
    ! Debitance
    real(rp) :: deb
    ! State variables for left crosssection
    real(rp) :: aL, debL, FrL, hL, kL, pL, rhL, uL, wL, zL, ddaL, dpdzL, dkdzL
    ! State variables for right crosssection
    real(rp) :: aR, debR, FrR, hR, kR, pR, rhR, uR, wR, zR, ddaR, dpdzR, dkdzR
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
    ! Coumpound channel state variables for left crosssection
    real(rp), dimension(3) :: acL, kcL, pcL, rhcL, wcL
    ! Coumpound channel state variables for right crosssection
    real(rp), dimension(3) :: acR, kcR, pcR, rhcR, wcR
#endif
    ! Residuals for convergence loop
    real(rp) :: errm1, errm2
    ! Discharge
    real(rp) :: q
    ! Water depths for convergence loop
    real(rp) :: h, hm1, hm2
    ! Critical water depth
    real(rp) :: hc
    ! Friction term
    real(rp) :: Sf
      
!======================================================================================================================!
!  External Functions
!======================================================================================================================!
    real(rp)  :: HtoS
      

#ifdef DEBUG_L1
    print *, "standard_step, segment:", iseg
#endif

    ! Retrieve state variables for last right cross-section (downstream cross-section)
    ics = mdl%msh%seg(iseg)%last_cs
    q = mdl%dof%q(ics)
    hR = mdl%dof%h(ics)
    zR = hR + mdl%msh%cs(ics)%bathy
    if (mdl%msh%strickler_type_code == strickler_type_average_conveyance .or. &
        mdl%msh%strickler_type_code == strickler_type_debord) then
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
        call areas_compound_channel(mdl%msh%cs(ics), hR, acR)
        aR = sum(acR)
        call perimeters_compound_channel(mdl%msh%cs(ics), hR, pcR)
        call widths_compound_channel(mdl%msh%cs(ics), hR, wcR)
        kcR(:) = mdl%msh%cs(ics)%strickler_params(:)
        do i = 1, 3
            if (pcR(i) > 1.e-12_rp) then
                rhcR(i) = acR(i) / pcR(i)
            else
                rhcR(i) = 0.d0
            end if
        end do
        aR = sum(acR)
        pR = sum(pcR)
        wR = sum(wcR)
        rhR = aR / pR
        debR = kcR(1) * acR(1) * rhcR(1)**d2p3 + &
               kcR(2) * acR(2) * rhcR(2)**d2p3 + &
               kcR(3) * acR(3) * rhcR(3)**d2p3
        ! print *, "OUT:H=", zR, ",A=", acR(2), acR(1)+acR(3), ",P=", pcR(2), pcR(1)+pcR(3)
#endif
    else
        call area(mdl%msh%cs(ics), hR, aR)
        call perimeter(mdl%msh%cs(ics), hR, pR)
        call width(mdl%msh%cs(ics), hR, wR)
        call strickler(mdl%msh%cs(ics), hR, kR)
        ! print *, "KEND", kR, mdl%msh%cs(ics)%strickler_params
        rhR = aR / pR
        debR = kR * aR * rhR**d2p3
        ! print *, "DEBEND", debR, aR, rhR
    end if
    uR = q / aR

      
    ! Iterate on cross-section from downstream to upstream
#ifdef DEBUG_SAVE_SEGMENT
    open(174, file="test_standard_step.csv")
    write(174, '(A)') "ics;dx;h;u;a;p;um_g;Sf;dxSf;dZdx;Q;QManning"
#endif
    do ics = mdl%msh%seg(iseg)%last_cs-1, mdl%msh%seg(iseg)%first_cs, -1

        !  Convergence loop
        hm1 = hR
        hm2 = hR + 0.1
        ! print *, "start:", hm1, hm2
        iter = 0
        do while(abs(hm1 - hm2) > eps .and. iter < itermax)
        
            iter = iter + 1
            mdl%dof%h(ics) = hm1
            zL = hm1 + mdl%msh%cs(ics)%bathy
            call update_level(mdl%msh%cs(ics), hm1)
          
            if (mdl%msh%strickler_type_code == strickler_type_average_conveyance .or. &
                mdl%msh%strickler_type_code == strickler_type_debord) then
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
                call areas_compound_channel(mdl%msh%cs(ics), hm1, acL)
                aL = sum(acL)
                call perimeters_compound_channel(mdl%msh%cs(ics), hm1, pcL)
                call widths_compound_channel(mdl%msh%cs(ics), hm1, wcL)
                kcL(:) = mdl%msh%cs(ics)%strickler_params(:)
                do i = 1, 3
                    if (pcL(i) > 1.e-12_rp) then
                        rhcL(i) = acL(i) / pcL(i)
                    else
                        rhcL(i) = 0.d0
                    end if
                end do
                aL = sum(acL)
                pL = sum(pcL)
                wL = sum(wcL)
                rhL = aL / pL
                debL = kcL(1) * acL(1) * rhcL(1)**d2p3 + &
                       kcL(2) * acL(2) * rhcL(2)**d2p3 + &
                       kcL(3) * acL(3) * rhcL(3)**d2p3
#endif
            else
                call area(mdl%msh%cs(ics), hm1, aL)
                call perimeter(mdl%msh%cs(ics), hm1, pL)
                call strickler(mdl%msh%cs(ics), hm1, kL)
                rhL = aL / pL
                debL = kL * aL * rhL**(d2p3)
            end if
            uL = q / aL
            
            deb = demi * (debL + debR)
            Sf = (q / deb)**2
#ifdef DEBUG_L2
            print *, "deb:", ics
            print *, "Detail:", debL, debR
            print *, "Detail+:", kL, rhL, aL, pL, hm1, abs(hm1-hm2), eps
            print *, "        ", kR, rhR, aR
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
            if (mdl%msh%strickler_type_code == strickler_type_average_conveyance .or. &
                mdl%msh%strickler_type_code == strickler_type_debord) then
                print *, "Detail:kcL", kcL(:)
                print *, "Detail:acL", acL(:)
                print *, "Detail:rhcL", rhcL(:)
                print *, "Detail:kcR", kcR(:)
                print *, "Detail:acR", acR(:)
                print *, "Detail:rhcR", rhcR(:)
                read(*, *)
            end if
#endif
            if (isnan(deb)) then
                print *, "deb is NaN:", ics, mdl%msh%cs(ics)%ibase
                print *, "Detail:", debL, debR
                print *, "Detail+:", kL, rhL, aL, pL, hm1
                print *, "        ", kR, rhR, aR
                read(*, *)
            end if
#endif
            call bisect_critical_depth(mdl%msh, mdl%dof, ics, mdl%gravity, q, hc)
            
            ! Compute elevation using energy equation
            dx = mdl%msh%cs(ics+1)%deltademi
            zL = zR + demi * (uR**2 - uL**2) / mdl%gravity + dx * Sf
            ! print *, demi * (uR**2 - uL**2) / mdl%gravity, dx * Sf, q, deb, dx
!             print *, " ** dx", dx
!             print *, " ** Hc", hc
!             read(*,*)
            
            h = zL - mdl%msh%cs(ics)%bathy
            
            ! Enforce minimal depth to critical depth
            call bisect_critical_depth(mdl%msh, mdl%dof, ics, mdl%gravity, q, hc)
            if (h < hc) h = hc
            if (h < mdl%heps) then
                h = mdl%heps
#ifndef CPP_ADJ
                mdl%warning_counters(1) = mdl%warning_counters(1) + 1
#endif
            end if
#ifdef DEBUG_L2
            print '(6(A,F12.6),A)', "L:d=", h, ",u=", uL, ",A=(", acL(2), ",", acL(1)+acL(3), &
                                  "),P=(", pcL(2), ",", pcL(1)+pcL(3), ")"
            print '(6(A,F12.6),A)', "R:d=", hR, ",u=", uR, ",A=(", acR(2), ",", acR(1)+acR(3), &
                                  "),P=(", pcR(2), ",", pcR(1)+pcR(3), ")"
            print *, " ** dZ", demi * (uR**2 - uL**2) / mdl%gravity + dx * Sf
            print *, " --DETAILS", uR, uL, Sf, dx
            print *, " ** z", zL
            print *, " ** DebL", debL, kcL(2) * acL(2) * rhcL(2)**d2p3, &
                                 kcL(1) * acL(1) * rhcL(1)**d2p3 + kcL(3) * acL(3) * rhcL(3)**d2p3
            print *, " ** DETAIL_L", kcL(2), acL(2), rhcL(2)
            print *, " ** DebR", debR, kcR(2) * acR(2) * rhcR(2)**d2p3, &
                                 kcR(1) * acR(1) * rhcR(1)**d2p3 + kcR(3) * acR(3) * rhcR(3)**d2p3
            print *, " ** DETAIL_R", kcR(2), acR(2), rhcR(2)
            print *, " ** h,hc", h, hc
            read(*,*)
#endif
            
            ! Update depth
            if (iter == 1) then
                errm1 = h - hm1
                hm1 = hm1 + 0.7 * errm1
            else
            
                errm2 = errm1
                errm1 = h - hm1
                if (abs(errm2 - errm1) < 1e-2) then
                    h = hm2 + demi * errm2
                else
                    h = hm2 - errm2 * (hm2 - hm1) / (errm2 - errm1)
                end if
            
                ! Enforce minimal depth to critical depth
                call bisect_critical_depth(mdl%msh, mdl%dof, ics, mdl%gravity, q, hc)
                if (h < hc) h = hc
                if (h < mdl%heps) then
                    h = mdl%heps
#ifndef CPP_ADJ
                    mdl%warning_counters(1) = mdl%warning_counters(1) + 1
#endif
                end if
                
                ! Prepare next step
                hm2 = hm1
                hm1 = h
                
                if (isnan(h)) then
                    mdl%status = 11
                    return
                    print *, " ** dZ", demi * (uR**2 - uL**2) / mdl%gravity + dx * Sf
                    print *, " --DETAILS", uR, uL, Sf, dx
                    print *, " ** h", hm2, hR
                    print *, " ** z", zL, zR
                    print *, " ** a", aL, aR
                    print *, " ** z", pL, pR
                    print *, " ** w", wL, wR
                    print *, " ** h,hc", h, hc
                    read(*,*)
                    print *, ics
                    print *, mdl%msh%cs(ics)%level_heights
                    print *, mdl%msh%cs(ics)%level_widths
                    read(*,*)
                end if
                
            end if
!             print *, " == ", h, hm1, hm2
                
        end do
            
        ! Compute final water depth
        h = demi * (hm1 + hm2)
        zL = mdl%msh%cs(ics)%bathy + h
        ! print *, ics, zL
        ! print *, "KCUR", kL, mdl%msh%cs(ics)%strickler_params
        ! print *, "DEBCUR", debL, aL, rhL
        
#ifdef DEBUG_SAVE_SEGMENT
        write(174, '(I4.4,11(A,E12.5))') ics-1, ";", dx, ";", h, ";", uL, ";", aL, ";", pL, &
                                        ";", demi * (uR**2 - uL**2) / mdl%gravity, ";", Sf, &
                                        ";", dx * Sf, ";", (zL - zR) /dx, ";", q, ";", debL * sqrt((zL - zR) /dx)
#endif

        mdl%dof%h(ics) = h
        mdl%dof%A(ics) = htoA(mdl%msh%cs(ics), mdl%dof%h(ics))
#ifdef DEBUG_L1
        print *, "INT", ics, mdl%dof%h(ics), mdl%dof%a(ics), mdl%dof%q(ics)
#endif

#ifndef CPP_ADJ
        ! Check for NaN value(s)
        if (isnan(mdl%dof%h(ics)) .or. isnan(mdl%dof%A(ics))) then
            mdl%status = 11
            return
        end if
#endif
!         print *, ics, mdl%dof%h(ics), mdl%msh%cs(ics)%bathy, kL, aL, uL
        
        ! Retrieve right state for next step
        hR = mdl%dof%h(ics)
        call update_level(mdl%msh%cs(ics), hR)
        zR = mdl%dof%h(ics) + mdl%msh%cs(ics)%bathy
        if (mdl%msh%strickler_type_code == strickler_type_average_conveyance .or. &
            mdl%msh%strickler_type_code == strickler_type_debord) then
#if defined(AVERAGE_CONVEYANCE) || defined(DEBORD_FORMULA)
            call areas_compound_channel(mdl%msh%cs(ics), hR, acR)
            aR = sum(acR)
            call perimeters_compound_channel(mdl%msh%cs(ics), hR, pcR)
            call widths_compound_channel(mdl%msh%cs(ics), hR, wcR)
            kcR(:) = mdl%msh%cs(ics)%strickler_params(:)
            do i = 1, 3
                if (pcR(i) > 1.e-12_rp) then
                    rhcR(i) = acR(i) / pcR(i)
                else
                    rhcR(i) = 0.d0
                end if
            end do
            aR = sum(acR)
            pR = sum(pcR)
            wR = sum(wcR)
            rhR = aR / pR
            debR = kcR(1) * acR(1) * rhcR(1)**d2p3 + &
                   kcR(2) * acR(2) * rhcR(2)**d2p3 + &
                   kcR(3) * acR(3) * rhcR(3)**d2p3
#endif
        else
            call area(mdl%msh%cs(ics), hR, aR)
            call perimeter(mdl%msh%cs(ics), hR, pR)
            call width(mdl%msh%cs(ics), hR, wR)
            call strickler(mdl%msh%cs(ics), hR, kR)
            rhR = aR / pR
            debR = kR * aR * rhR**d2p3
        end if
        uR = q / aR
!         print *, "STANDARD_STEP(0):", ics, zR, aR, kR, mdl%msh%cs(ics)%strickler_params
!         if (ics ==  mdl%msh%seg(iseg)%first_cs) read(*,*)
        
    end do

#ifdef DEBUG_SAVE_SEGMENT
    close(174)
#endif

end subroutine
