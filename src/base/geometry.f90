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
!> \file geometry.f90
!> \brief This file includes all the geometry routines

!> Compute polynomials approximations for a cross-section
subroutine compute_levels_poly(cs)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param cs Cross-section
    type(Crosssection), intent(inout) :: cs
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of level k
    integer(ip) :: ik
    ! Degree of polynomial
    integer(ip) :: degree
    ! Delta area between k-1 and k
    real(rp) :: delta_area
    ! Delta perimeter between k-1 and k
    real(rp) :: delta_perim
    ! Height at level k
    real(rp) :: Zk
    ! Height at level k-1
    real(rp) :: Zkm1
    ! Width at level k
    real(rp) :: Wk
    ! Width at level k-1
    real(rp) :: Wkm1
  
    degree = size(cs%poly, 1) - 1
    select case(degree)
    
        case(1)
        
            cs%poly(1, 1) = 0.0
            cs%poly(2, 1) = cs%level_widths(1)
            
            do ik = 2, cs%nlevels
            
                Zk = cs%level_heights(ik)
                Zkm1 = cs%level_heights(ik-1)
                Wk = cs%level_widths(ik)
                Wkm1 = cs%level_widths(ik-1)
                cs%poly(1, ik) = (Wk - Wkm1) / (Zk - Zkm1)
                cs%poly(2, ik) = Wk - cs%poly(1, ik) * Zk
            end do
            
!         case(3)
!         
!             call f90wrap_abort("Computing polynoms coefficients for degree 3 is not implemented yet")
!             
!         case default
!         
!             call f90wrap_abort("degree must be 1 or 3 to compute polynoms coefficients")
            
    end select
      
end subroutine


subroutine compute_levels_cum(cs)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Index of quadrature point
    integer(ip) :: iq
    ! Degree of polynomial
    integer(ip) :: degree
    ! Delta area between k-1 and k
    real(rp) :: delta_area
    ! Delta pa between k-1 and k
    real(rp) :: delta_pa
    ! Delta perimeter between k-1 and k
    real(rp) :: delta_perim
    ! Gradient of width w.r.t. height
    real(rp) :: dwdz
    ! Height
    real(rp) :: z
    ! Height at level k
    real(rp) :: Zk
    ! Height at level k-1
    real(rp) :: Zkm1
    ! Width at level k
    real(rp) :: Wk
    ! Width at level k-1
    real(rp) :: Wkm1
    !=PARAMETERS=======================================================================================================!
    ! Quadrature points (for polynomials up to degree 5)
    real(rp), dimension(2), parameter :: ksiq = (/-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)/)
    ! Quadrature weights (for polynomials up to degree 5)
    real(rp), dimension(2), parameter :: wq = (/1.0, 1.0/)
  
    ! Compute cumulated flow area
    cs%area_cum(0) = 0.0
    cs%area_cum(1) = (cs%level_heights(1) - cs%bathy) * cs%level_widths(1)
!     print *, "levels_cum"
!     print *, "aera_cum(1)=", cs%area_cum(1)
    do ik = 2, cs%nlevels
    
!         print *, "ik=", ik
        delta_area = 0.0
        do ic = 1, size(cs%poly, 1)
            degree = size(cs%poly, 1) - ic
            Zk = cs%level_heights(ik)
            Zkm1 = cs%level_heights(ik-1)
            delta_area = delta_area + cs%poly(ic, ik) / (degree + 1) * (Zk**(degree+1) - Zkm1**(degree+1))
!             print *, "delta_area+=", cs%poly(ic, ik) / (degree + 1) * Zk**(degree+1), "-", &
!                                      cs%poly(ic, ik) / (degree + 1) * Zkm1**(degree+1)
        end do
        
        cs%area_cum(ik) = cs%area_cum(ik-1) + delta_area
!         print *, "=>=", cs%area_cum(ik), cs%area_cum(ik-1), "+", delta_area
    end do 

    ! Compute cumulated wetted perimeter
    cs%perim_cum(0) = cs%level_widths(1)
    cs%perim_cum(1) = 2 * (cs%level_heights(1) - cs%bathy) + cs%level_widths(1)
    degree = size(cs%poly, 1) - 1
    if (degree == 1) then
        do ik = 2, cs%nlevels
            Zk = cs%level_heights(ik)
            Zkm1 = cs%level_heights(ik-1)
            Wk = cs%level_widths(ik)
            Wkm1 = cs%level_widths(ik-1)
            delta_perim = 2.0 * sqrt((Zk - Zkm1)**2 + (0.5 * (Wk - Wkm1))**2)
            cs%perim_cum(ik) = cs%perim_cum(ik-1) + delta_perim
        end do
    else if (degree == 2) then
        do ik = 2, cs%nlevels
            Zk = cs%level_heights(ik)
            Zkm1 = cs%level_heights(ik-1)
            
            ! Compute integral using Gauss-Legendre quadrature
            delta_perim = 0.0
            do iq = 1, 2
                z = 0.5 * (Zk - Zkm1) * ksiq(iq) + 0.5 * (Zkm1 + Zk)
                dwdz = 2.0 * cs%poly(1, ik) * z + cs%poly(2, ik)
                delta_perim = delta_perim + wq(iq) * sqrt(1.0 + dwdz**2)
            end do
            delta_perim = delta_perim * 0.5 * (Zk - Zkm1)
            
            cs%perim_cum(ik) = cs%perim_cum(ik-1) + delta_perim
            
        end do
    else if (degree == 3) then
        do ik = 2, cs%nlevels
            Zk = cs%level_heights(ik)
            Zkm1 = cs%level_heights(ik-1)
            
            ! Compute integral using Gauss-Legendre quadrature
            delta_perim = 0.0
            do iq = 1, 2
                z = 0.5 * (Zk - Zkm1) * ksiq(iq) + 0.5 * (Zkm1 + Zk)
                dwdz = 3.0 * cs%poly(1, ik) * z**2 + 2.0 * cs%poly(2, ik) * z + cs%poly(3, ik)
                delta_perim = delta_perim + wq(iq) * sqrt(1.0 + dwdz**2)
            end do
            delta_perim = delta_perim * 0.5 * (Zk - Zkm1)
            
            cs%perim_cum(ik) = cs%perim_cum(ik-1) + delta_perim
            
        end do
!     else
!         call f90wrap_abort("Degree must be 1 or 3 to compute cumulated wetted perimeter")
    end if
  
    ! Compute cumulated pa
    do ik = 1, cs%nlevels
    
        delta_pa = 0.0
        do ic = 1, size(cs%poly, 1)
            degree = size(cs%poly, 1) - ic
            Zk = cs%level_heights(ik)
            delta_pa = delta_pa + cs%poly(ic, ik) / (degree + 1) * (Zk**(degree+1) - cs%bathy**(degree+1))
        end do
        
        if (ik == 1) then
            cs%pa_cum(ik) = delta_pa
        else
            cs%pa_cum(ik) = cs%pa_cum(ik-1) + delta_pa
        end if
    end do 
      
end subroutine
       
    
!> Update the level of a cross-section given a depth
subroutine update_level(cs, h)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h

    if (h < 0.0) then
        cs%level = 0
    else if (h <= cs%level_heights(1) - cs%bathy) then
        cs%level = 1
#ifdef FLOODPLAIN_MODEL
    else if (h >= cs%level_heights(cs%nlevels) - cs%bathy) then
        cs%level = cs%nlevels + 1
#else
    else if (h > cs%level_heights(cs%nlevels) - cs%bathy) then
        cs%level = cs%nlevels
#endif
    else
        if (cs%level < 1) cs%level = 1
        if (cs%level > cs%nlevels) cs%level = cs%nlevels
        if (cs%level > 1) then
            do while(h + cs%bathy < cs%level_heights(cs%level))
                cs%level = cs%level - 1
                if (cs%level < 1) exit
            end do
        end if
        if (cs%level < cs%nlevels) then
            do while(h + cs%bathy > cs%level_heights(cs%level))
                cs%level = cs%level + 1
                if (cs%level >= cs%nlevels) exit
            end do
        end if
    end if
  
end subroutine
       
    
!> Compute the height of a cross-section given a depth
subroutine update_level_from_area(cs, area)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: area

    if (area < 0.0) then
        cs%level = 0
    else if (area < cs%area_cum(1)) then
        cs%level = 1
#ifdef FLOODPLAIN_MODEL
    else if (area > cs%area_cum(cs%nlevels)) then
        cs%level = cs%nlevels + 1
#else
    else if (area > cs%area_cum(cs%nlevels)) then
        cs%level = cs%nlevels
#endif
    else
        if (cs%level < 1) cs%level = 1
        if (cs%level > 1) then
            do while(area < cs%area_cum(cs%level))
                cs%level = cs%level - 1
            end do
        end if
        if (cs%level < cs%nlevels) then
            do while(area > cs%area_cum(cs%level))
                cs%level = cs%level + 1
            end do
        end if
    end if
  
end subroutine
       
    
!> Compute the area of a cross-section given a depth
function htoA(cs, h) result(A)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h
    real(rp) :: A
    
    call update_level(cs, h)
    call area(cs, h, A)

end function

       
    
!> Compute the area of a cross-section given a depth without updating the level
function htoA_noupdate(cs, h) result(A)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h
    real(rp) :: A
    !=EXTERNAL FUNCTIONS===============================================================================================!
    real(rp) :: htoA
    !=LOCAL VARIABLES==================================================================================================!
    integer(ip) :: level

    ! Store current level
    level = cs%level
    
    ! Call htoA function
    A = htoA(cs, h)

    ! Restore current level if noupdate is true
    cs%level = level

end function htoA_noupdate

       
    
!> Compute the wetted perimeter of a cross-section given a depth
function htoP(cs, h) result(P)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h
    real(rp) :: P

    call update_level(cs, h)
    call perimeter(cs, h, P)
  
end function

       
    
!> Compute the wetted perimeter of a cross-section given a depth without updating the level
function htoP_noupdate(cs, h) result(P)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h
    real(rp) :: P
    !=EXTERNAL FUNCTIONS===============================================================================================!
    real(rp) :: htoP
    !=LOCAL VARIABLES==================================================================================================!
    !> Index of level
    integer(ip) :: level

    ! Store current level
    level = cs%level
    
    ! Call htoP function
    P = htoP(cs, h)

    ! Restore current level if noupdate is true
    cs%level = level

end function htoP_noupdate
       
    
!> Compute the top width of a cross-section given a depth
function htoW(cs, h) result(W)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h
    real(rp) :: W
    
    call update_level(cs, h)
    call width(cs, h, W)
  
end function

       
    
!> Compute the area of a cross-section given a depth without updating the level
function htoW_noupdate(cs, h) result(W)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: h
    real(rp) :: W
    !=EXTERNAL FUNCTIONS===============================================================================================!
    real(rp) :: htoW
    !=LOCAL VARIABLES==================================================================================================!
    !> Index of level
    integer(ip) :: level

    ! Store current level
    level = cs%level
    
    ! Call htoW function
    W = htoW(cs, h)

    ! Restore current level if noupdate is true
    cs%level = level

end function htoW_noupdate

    
!> @todo
function Atoh(cs, A) result(h)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: A
    real(rp) :: h

    call update_level_from_area(cs, A)
    call depth_from_area(cs, A, h)
  
end function

       
    
!> Compute the area of a cross-section given a depth without updating the level
function Atoh_noupdate(cs, A) result(h)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(inout) :: cs
    real(rp), intent(in) :: A
    real(rp) :: h
    !=EXTERNAL FUNCTIONS===============================================================================================!
    real(rp) :: Atoh
    !=LOCAL VARIABLES==================================================================================================!
    !> Index of level
    integer(ip) :: level

    ! Store current level
    level = cs%level
    
    ! Call Atoh function
    h = Atoh(cs, A)

    ! Restore current level if noupdate is true
    cs%level = level

end function Atoh_noupdate

    
!> Compute the wetted area of a cross-section given a depth
subroutine width(cs, h, value)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Degree of coefficient
    integer(ip) :: degree
    ! Free surface elevation
    real(rp) :: Z
    ! Free surface elevation at level k-1
    real(rp) :: Zkm1

    ik = cs%level
    if (ik == 0) then
        value = cs%level_widths(1)
#ifdef FLOODPLAIN_MODEL
    else if (ik > cs%nlevels) then
        Z = h + cs%bathy
!         Zk = cs%level_heights(cs%nlevels)
        Zkm1 = cs%level_heights(cs%nlevels)
        value = cs%level_widths(cs%nlevels) + (Z - Zkm1) * cs%alphaFP
#endif
    else
        value = 0.0
        Z = h + cs%bathy
        do ic = 1, size(cs%poly, 1)
            degree = size(cs%poly, 1) - ic
            value = value + cs%poly(ic, ik) * Z**degree
        end do  
    end if
  
end subroutine

    
!> Compute the wetted area of a cross-section given a depth
subroutine widths_compound_channel(cs, h, values)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    real(rp), dimension(3), intent(out) :: values
    !=LOCAL VARIABLES==================================================================================================!
    ! Total width
    real(rp) :: total_width
    
    call width(cs, h, total_width)
    
    if (cs%ob_levels(1) > 0 .and. cs%ob_levels(2) > 0) then
        values(1) = max(0.0_rp, 0.5 * (total_width - cs%level_widths(cs%ob_levels(1))))
        values(3) = max(0.0_rp, 0.5 * (total_width - cs%level_widths(cs%ob_levels(2))))
    else if (cs%ob_levels(1) > 0) then
        values(1) = max(0.0_rp, 0.5 * (total_width - cs%level_widths(cs%ob_levels(1))))
        values(3) = 0.0
    else if (cs%ob_levels(2) > 0) then
        values(1) = 0.0
        values(3) = max(0.0_rp, 0.5 * (total_width - cs%level_widths(cs%ob_levels(2))))
    else
        values(1) = 0.0
        values(3) = 0.0
    end if
    values(2) = total_width - values(1) - values(3)
  
end subroutine

    
!> Compute the wetted area of a cross-section given a depth
subroutine depth_from_area(cs, area, h)
    ! TODO
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: area
    real(rp), intent(out) :: h
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Degree of coefficient
    integer(ip) :: degree
    ! Variables for roots evaluation for a polynomial of degree 2
    real(rp) :: a, b, c, delta, root1, root2
    ! Free surface elevation
    real(rp) :: Z
    ! Free surface elevation at level k-1
    real(rp) :: Zkm1

    ik = cs%level
    if (ik == 0) then
        h = 0.0
#ifdef FLOODPLAIN_MODEL
    else if (ik > cs%nlevels) then
        a = cs%alphaFP
        b = cs%level_widths(cs%nlevels)
        c = cs%area_cum(cs%nlevels) - area
        delta = b**2 - 4.0 * a * c
        root1 = (-b + sqrt(delta)) / (2.0 * a)
        Zkm1 = cs%level_heights(ik)
        h = root1 + Zkm1 - cs%bathy
#endif
    else
        degree = size(cs%poly, 1) - 1
!         print *, "depth_from_area:area=", area, ik
        if (degree == 1) then
            if (ik == 1) then
                Zkm1 = cs%bathy
            else
                Zkm1 = cs%level_heights(ik-1)
            end if
            if (abs(cs%poly(1, ik)) < 1e-12) then
!                 Zkm1 = cs%level_heights(ik-1)
                Z = (area - cs%area_cum(ik-1) + cs%poly(2, ik) * Zkm1) / cs%poly(2, ik)
                h = Z - cs%bathy
                return
            else
!                 if (ik == 1) then
!                     Zkm1 = cs%bathy
!                 else
!                     Zkm1 = cs%level_heights(ik-1)
!                 end if
!                 Zk = cs%level_heights(ik-1)
                a = 0.5 * cs%poly(1, ik)
                b = cs%poly(2, ik)
                c = cs%area_cum(ik-1) - area - 0.5 * cs%poly(1, ik) * Zkm1**2 - cs%poly(2, ik) * Zkm1
                delta = b**2 - 4.0 * a * c
                root1 = (-b + sqrt(delta)) / (2.0 * a)
                root2 = (-b - sqrt(delta)) / (2.0 * a)
                h = root1 - cs%bathy
!                 if (root2 > Zk) then
!                   h = root2 - cs%bathy
!                 else
!                   h = root1 - cs%bathy
!                 end if
!                 print *, "--roots=", root1, root2, Zk
            end if
        end if
    end if
  
end subroutine

    
!> Compute the wetted area of a cross-section given a depth
subroutine area(cs, h, value)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Degree of coefficient
    integer(ip) :: degree
    ! Top width
    real(rp) :: W
    ! Free surface elevation
    real(rp) :: Z
    ! Free surface elevation at level k-1
    real(rp) :: Zkm1

    ik = cs%level
    if (ik == 0) then
        value = 0.0
#ifdef FLOODPLAIN_MODEL
    else if (ik > cs%nlevels) then
        Z = h + cs%bathy
        Zkm1 = cs%level_heights(cs%nlevels)
        W = cs%level_widths(cs%nlevels) + 2.0 * cs%alphaFP * (Z - Zkm1)
        value = cs%area_cum(cs%nlevels) + (Z - Zkm1) * 0.5 * (W + cs%level_widths(cs%nlevels))
! #ifndef CPP_ADJ
!         if (value < 0.0) then
!             print *, "area < 0 (FLOODPLAIN_MODEL)"
!             print *, "Detail:", ik, cs%nlevels
!             print *, "Detail+:", W, cs%alphaFP, Z, Zkm1
!             read(*, *)
!         end if
! #endif
!         print *, "[FP] area=", value, cs%area_cum(cs%nlevels)
!         print *, "--Z:", Z
!         print *, "--Zk:", Zk
!         print *, "--W:", W
!         print *, "--Wk:", cs%level_widths(cs%nlevels)
!         read(*,*)
#else
    else if (ik > cs%nlevels) then
        Z = h + cs%bathy
        Zkm1 = cs%level_heights(cs%nlevels)
        value = cs%area_cum(cs%nlevels) + (Z - Zkm1) * cs%level_widths(cs%nlevels)
#endif
    else
        value = 0.0
        Z = h + cs%bathy
        if (ik == 1) then
            Zkm1 = cs%bathy
        else
            Zkm1 = cs%level_heights(ik-1)
        end if
!         value = h * cs%level_widths(1)
        do ic = 1, size(cs%poly, 1)
            degree = size(cs%poly, 1) - ic
!             Zkm1 = cs%level_heights(ik-1)
!             value = value + cs%poly(ic, ik) / (degree + 1) * (Z**(degree+1) - Zk**(degree+1))
            value = value + cs%poly(ic, ik) / (degree + 1) * (Z**(degree+1) - Zkm1**(degree+1))
        end do  
        if (ik > 1) value = cs%area_cum(ik-1) + value
! #ifndef CPP_ADJ
!         if (value < 1e-12 .and. h > 1e-12) then
!             print *, "area < 0 (POLY)", value
!             print *, "Detail:", ik, cs%nlevels
!             print *, "Detail+:", cs%poly(:, ik)
!             print *, "bathy:", cs%bathy
!             print *, "level_heights:", cs%level_heights
!             print *, "level_widths:", cs%level_widths
!             Z = h + cs%bathy
!             print *, "Serial:", Z
!             value = 0.0
!             do ic = 1, size(cs%poly, 1)
!                 degree = size(cs%poly, 1) - ic
! !                 Zk = cs%level_heights(ik)
! !                 value = value + cs%poly(ic, ik) / (degree + 1) * (Z**(degree+1) - Zk**(degree+1))
!                 value = value + cs%poly(ic, ik) / (degree + 1) * (Z**(degree+1) - Zkm1**(degree+1))
!                 print *, "value+=:", cs%poly(ic, ik), "/", (degree + 1), "*(", Z, "**", (degree+1), "-", &
!                                      Zkm1, "**", (degree+1), ")"
!             end do  
!             read(*, *)
!         end if
! #endif
    end if
  
end subroutine

    
!> Compute the wetted areas of the three parts (left overbank, channel, right overbank) of a cross-section, given a
!> depth
subroutine areas_compound_channel(cs, h, values)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    real(rp), dimension(3), intent(out) :: values
    !=LOCAL VARIABLES==================================================================================================!
    ! Delta depth
    real(rp) :: dh
    ! Rectangle area
    real(rp) :: rect_area
    ! Total wetted area
    real(rp) :: total_area
    
    call area(cs, h, total_area)
    
!     print *, "areas_compound_channel:h=", h
    
    if (cs%ob_levels(1) > 0 .and. cs%ob_levels(2) > 0) then
        dh = h - (cs%level_heights(cs%ob_levels(1)) - cs%bathy)
        rect_area = max(0.0_rp, cs%level_widths(cs%ob_levels(1)) * dh)
        values(1) = max(0.0_rp, 0.5 * (total_area - cs%area_cum(cs%ob_levels(1)) - rect_area))
        dh = h - (cs%level_heights(cs%ob_levels(2)) - cs%bathy)
        rect_area = max(0.0_rp, cs%level_widths(cs%ob_levels(2)) * dh)
        values(3) = max(0.0_rp, 0.5 * (total_area - cs%area_cum(cs%ob_levels(2)) - rect_area))
    else if (cs%ob_levels(1) > 0) then
        dh = h - (cs%level_heights(cs%ob_levels(1)) - cs%bathy)
        rect_area = max(0.0_rp, cs%level_widths(cs%ob_levels(1)) * dh)
        values(1) = max(0.0_rp, 0.5 * (total_area - cs%area_cum(cs%ob_levels(1)) - rect_area))
        values(3) = 0.0
    else if (cs%ob_levels(2) > 0) then
        values(1) = 0.0
        dh = h - (cs%level_heights(cs%ob_levels(2)) - cs%bathy)
        rect_area = max(0.0_rp, cs%level_widths(cs%ob_levels(2)) * dh)
        values(3) = max(0.0_rp, 0.5 * (total_area - cs%area_cum(cs%ob_levels(2)) - rect_area))
    else
        values(1) = 0.0
        values(3) = 0.0
    end if
    values(2) = total_area - values(1) - values(3)
!     print *, "-LOB:A=", values(1)
!     print *, "-MC :A=", values(2)
!     print *, "-ROB:A=", values(3)
  
end subroutine


!> Compute the wetted perimeter of a cross-section given a depth
subroutine perimeter(cs, h, value)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Index of quadrature point
    integer(ip) :: iq
    ! Power coefficient
    integer(ip) :: power
    ! Degree of polynomial
    integer(ip) :: degree
    ! Gradient of W w.r.t. Z
    real(rp) :: dwdz
    ! Free surface width
    real(rp) :: W
    ! Free surface width at level k
    real(rp) :: Wk
    ! Free surface elevation
    real(rp) :: Z
    ! Free surface elevation at level k-1
    real(rp) :: Zkm1
    ! Free surface elevation
    real(rp) :: Zq
    !=PARAMETERS=======================================================================================================!
    ! Quadrature points (for polynomials up to degree 5)
    real(rp), dimension(2), parameter :: ksiq = (/-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)/)
    ! Quadrature weights (for polynomials up to degree 5)
    real(rp), dimension(2), parameter :: wq = (/1.0, 1.0/)

    ik = cs%level
    if (ik == 0) then
        value = cs%level_widths(1)
    else if (ik == 1) then
        value = 2.0 * h + cs%level_widths(1)
#ifdef FLOODPLAIN_MODEL
    else if (ik > cs%nlevels) then
        Z = h + cs%bathy
        Zkm1 = cs%level_heights(cs%nlevels)
        W = cs%level_widths(cs%nlevels) + 2.0 * cs%alphaFP * (Z - Zkm1)
        value = cs%perim_cum(cs%nlevels) + 2.0 * (Z - Zkm1) * sqrt(1.0 + cs%alphaFP**2)
#endif
    else
        value = 0.0
        Z = h + cs%bathy
        if (ik == 1) then
            Zkm1 = cs%bathy
        else
            Zkm1 = cs%level_heights(ik-1)
        end if
!         Zk = cs%level_heights(ik)
        degree = size(cs%poly, 1) - 1
        if (degree == 1) then
            value = 2.0 * sqrt(1.0 + (0.5 * cs%poly(1, ik))**2) * (Z - Zkm1)
!             print *, "perimeter(deg1, ik>0)::", h, value
        else if (degree < 5) then
          
            ! Compute int(Zk,Z)[sqrt(1 + dwdz)dz] using Gauss-Legendre quadrature
            ! n=2 suitable for polynomials of degree 2n-1 <= 5
            do iq = 1, 2
              Zq = 0.5 * (Z - Zkm1) * ksiq(iq) + 0.5 * (Z - Zkm1)
              dwdz = 0.0
              do ic = 1, size(cs%poly, 1)-1
                power = size(cs%poly, 1) - ic
                dwdz = dwdz + real(power, rp) * cs%poly(ic, ik) * Zq**(power-1)
              end do
              value = value + wq(iq) * sqrt(1.0 + dwdz**2)
            end do
            value = value * 0.5 * (Z - Zkm1)
            
        end if
        
        ! Perimeter is equal to cumulated perimeter on [1, ik-1] + int(Z,Zk-1)[sqrt(1 + dwdz)dz]
        value = cs%perim_cum(ik-1) + value
        
    end if
  
end subroutine


!> Compute the wetted perimeters of the three parts (left overbank, channel, right overbank) of a cross-section, given a
!> depth
subroutine perimeters_compound_channel(cs, h, values)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    !> Perimeters (left overbank, channel, right overbank)
    real(rp), dimension(3), intent(out) :: values
    !=LOCAL VARIABLES==================================================================================================!
    ! Total wetted perimeter
    real(rp) :: total_perimeter
    
    call perimeter(cs, h, total_perimeter)
    
    if (cs%ob_levels(1) > 0 .and. cs%ob_levels(2) > 0) then
        values(1) = max(0.0_rp, 0.5 * (total_perimeter - cs%perim_cum(cs%ob_levels(1))))
        values(3) = max(0.0_rp, 0.5 * (total_perimeter - cs%perim_cum(cs%ob_levels(2))))
    else if (cs%ob_levels(1) > 0) then
        values(1) = max(0.0_rp, total_perimeter - cs%perim_cum(cs%ob_levels(1)))
        values(3) = 0.0
    else if (cs%ob_levels(2) > 0) then
        values(1) = 0.0
        values(3) = max(0.0_rp, total_perimeter - cs%perim_cum(cs%ob_levels(2)))
    else
        values(1) = 0.0
        values(3) = 0.0
    end if
    values(2) = total_perimeter - values(1) - values(3)
  
end subroutine

    
!> Compute the Strickler coefficient of a cross-section given a depth
subroutine strickler(cs, h, value)
    use m_common
    use m_mesh, only: Crosssection
    use m_mesh, only: strickler_type_constant
    use m_mesh, only: strickler_type_powerlaw_h
#ifdef STRICKLER_EINSTEIN
    use m_mesh, only: strickler_type_einstein
#endif
#ifdef AVERAGE_CONVEYANCE
    use m_mesh, only: strickler_type_average_conveyance
#endif
#ifdef DEBORD_FORMULA
    use m_mesh, only: strickler_type_debord
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Cross-section
    type(Crosssection), intent(in) :: cs
    !> Water depth
    real(rp), intent(in) :: h
    !> Computed Strickler
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
#ifdef STRICKLER_EINSTEIN
    ! Perimeters of compound channel
    real(rp), dimension(3) :: P
    ! Computation variables
    real(rp) :: A, B
#endif
        
    if (cs%strickler_type_code == strickler_type_constant) then
        value = cs%strickler_params(1)
    else if (cs%strickler_type_code == strickler_type_powerlaw_h) then
        value = cs%strickler_params(1) * h**cs%strickler_params(2)
#ifdef STRICKLER_EINSTEIN
    else if (cs%strickler_type_code == strickler_type_einstein) then
    
        ! Compute partial wetted perimeters
        call perimeters_compound_channel(cs, h, P)
!         print *, "P_compound_channel:", P
!         print *, "K_compound_channel:", cs%strickler_params(:)
        
        A = P(1) + P(2) + P(3)
        B = P(1) / cs%strickler_params(1)**d3p2 + P(2) / cs%strickler_params(2)**d3p2 + &
            P(3) / cs%strickler_params(3)**d3p2
        value = (A / B)**d2p3
#endif
#ifdef AVERAGE_CONVEYANCE
    else if (cs%strickler_type_code == strickler_type_average_conveyance) then
#ifndef CPP_ADJ
        call f90wrap_abort("Subroutine 'strickler' must not be called with 'average_conveyance' model")
#endif
#endif
#ifdef DEBORD_FORMULA
    else if (cs%strickler_type_code == strickler_type_debord) then
#ifndef CPP_ADJ
        call f90wrap_abort("Subroutine 'strickler' must not be called with 'Debord' model")
#endif
#endif
    end if
  
end subroutine

    
!> Compute the pressure of a cross-section given a depth
subroutine pressure(cs, h, gravity, value)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Crosssection), intent(in) :: cs
    real(rp), intent(in) :: h
    real(rp), intent(in) :: gravity
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Degree of coefficient
    integer(ip) :: degree
    ! Multipliers
    integer(ip) :: mult1, mult2
    ! Free surface elevation
    real(rp) :: Z
!     ! Free surface elevation at level k-1
!     real(rp) :: Zkm1

    ik = cs%level
    if (ik == 0) then
        value = 0.0
!         value = gravity * 0.5 * h**2 * cs%level_widths(1)
    else
        value = 0.0
        Z = h + cs%bathy
        mult1 =  1.0 / (degree + 1) * (degree + 2)
        mult2 =  1.0 / (degree + 1)
!         if (ik == 1) then
!             Zkm1 = cs%bathy
!         else
!             Zkm1 = cs%level_heights(ik-1)
!         end if
        do ic = 1, size(cs%poly, 1)
            degree = size(cs%poly, 1) - ic
!             Zk = cs%level_heights(ik)
            value = value + cs%poly(ic, ik) * mult1 * (Z**(degree+2) - cs%bathy**(degree+2))
            value = value - cs%poly(ic, ik) * mult2 * cs%bathy**(degree+1) * h
        end do  
        value = gravity * (cs%pa_cum(ik) * h + value)
    end if
  
end subroutine

    
!> Compute the Sg term of a cross-section given a depth
subroutine Sg(cs, csm1, csp1, h, gravity, value)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Current cross-section
    type(Crosssection), intent(in) :: cs
    !> Upstream cross-section (i-1)
    type(Crosssection), intent(in) :: csm1
    !> Downstream cross-section (i+1)
    type(Crosssection), intent(in) :: csp1
    real(rp), intent(in) :: h
    real(rp), intent(in) :: gravity
    real(rp), intent(out) :: value
    
    ! TODO
!     ! Index of polynomial coefficient
!     integer(ip) :: ic
!     ! Index of level k
!     integer(ip) :: ik
!     ! Degree of coefficient
!     integer(ip) :: degree
!     ! Multipliers
!     integer(ip) :: mult1, mult2
!     ! Free surface elevation
!     real(rp) :: Z
!     ! Free surface elevation at level k
!     real(rp) :: Zk
!     
!     ! Store levels for section i-1 and i+1
!     level_im1 = csm1%level
!     level_ip1 = csp1%level
!     
!     ! Update levels for cross-sections i-1 and i+1 using hi
!     call update_level(csm1, h)
!     call update_level(csp1, h)
! 
!     ! Compute pressure at i-1 and i+1
!     call pressure(csm1, h, P_im1)
!     call pressure(csm1, h, P_ip1)
!     
!     
!     if (ik == 0) then
!         pressure = gravity * 0.5 * h**2 * cs%level_widths(1)
!     else
!         pressure = 0.0
!         Z = h + cs%bathy
!         mult1 =  1.0 / (degree + 1) * (degree + 2)
!         mult2 =  1.0 / (degree + 1)
!         do ic = 1, size(cs%poly, 1)
!             degree = size(cs%poly, 1) - ic
!             Zk = cs%level_heights(ik)
!             pressure = pressure + cs%poly(ic, ik) * mult1 * (Z**(degree+2) - cs%bathy**(degree+2))
!             pressure = pressure - cs%poly(ic, ik) * mult2 * cs%bathy**(degree+1) * h
!         end do  
!         pressure = gravity * (cs%pa_cum(ik) * h + pressure)
!     end if
  
end subroutine

    
!> Compute the gradient of P w.r.t. Z of a cross-section given the depth
subroutine dPdZ(cs, h, value)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Cross-section
    type(Crosssection), intent(in) :: cs
    !> Water depth
    real(rp), intent(in) :: h
    !> Computed Strickler
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of polynomial coefficient
    integer(ip) :: ic
    ! Index of level k
    integer(ip) :: ik
    ! Degree of polynomial
    integer(ip) :: degree
    ! Gradient of W w.r.t. Z
    real(rp) :: dwdz
    ! Power coefficient
    integer(ip) :: power
    ! Free surface elevation
    real(rp) :: Z

    ik = cs%level
    if (ik == 0) then
        value = 1.0
#ifdef FLOODPLAIN_MODEL
    else if (ik > cs%nlevels) then
        value = 2.0 * sqrt(1.0 + cs%alphaFP**2)
#endif
    else
        degree = size(cs%poly, 1) - 1
        if (degree == 1) then
            value = 2.0 * sqrt(1.0 + (0.5 * cs%poly(1, ik))**2)
        else
            Z = h + cs%bathy
            dwdz = 0.0
            do ic = 1, size(cs%poly, 1)-1
              power = size(cs%poly, 1) - ic
              dwdz = dwdz + real(power, rp) * cs%poly(ic, ik) * Z**(power-1)
            end do
            value = sqrt(1.0 + dwdz**2)
            
        end if
        
    end if
    
end subroutine

    
!> Compute the gradient of three parts perimeters of compound channel w.r.t. Z of a cross-section, given the depth
subroutine dPdZ_compound_channel(cs, h, values)
    use m_common
    use m_mesh, only: Crosssection
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Cross-section
    type(Crosssection) :: cs
    !> Water depth
    real(rp), intent(in) :: h
    !> Computed dPdZ (left overbank, channel, right overbank)
    real(rp), dimension(3), intent(out) :: values
    !=LOCAL VARIABLES==================================================================================================!
    ! Level
    integer(ip) :: level
    real(rp), parameter :: dh = 0.1
    ! Total dPdZ
    real(rp) :: total_dPdZ
    ! Total wetted perimeter
    real(rp) :: total_perimeter
    ! Wetted perimeters for gradient computation
    real(rp), dimension(3) :: Pp, Pm
    
    level = cs%level
    call update_level(cs, h+dh)
    call perimeters_compound_channel(cs, h+dh, Pp)
    call update_level(cs, h-dh)
    call perimeters_compound_channel(cs, h-dh, Pm)
    cs%level = level
    values(:) = (Pp(:) - Pm(:)) / (2 * dh)
!     
!     if (cs%ob_levels(1) > 0 .and. cs%ob_levels(2) > 0) then
!         if (cs%perim_cum(cs%ob_levels(1)) < total_perimeter .and. &
!             cs%perim_cum(cs%ob_levels(2)) < total_perimeter) then
!             values(1) = 0.5 * total_dPdZ
!             values(2) = 0.0
!             values(3) = 0.5 * total_dPdZ
!         else if (cs%perim_cum(cs%ob_levels(1)) < total_perimeter) then
!             values(1) = total_dPdZ
!             values(2) = 0.0
!             values(3) = 0.0
!         else if (cs%perim_cum(cs%ob_levels(2)) < total_perimeter) then
!             values(1) = 0.0
!             values(2) = 0.0
!             values(3) = total_dPdZ
!         else
!             values(1) = 0.0
!             values(2) = total_dPdZ
!             values(3) = 0.0
!         end if
!     else if (cs%ob_levels(1) > 0) then
!         if (cs%perim_cum(cs%ob_levels(1)) < total_perimeter) then
!             values(1) = total_dPdZ
!             values(2) = 0.0
!             values(3) = 0.0
!         else
!             values(1) = 0.0
!             values(2) = total_dPdZ
!             values(3) = 0.0
!         end if
!     else if (cs%ob_levels(2) > 0) then
!         if (cs%perim_cum(cs%ob_levels(2)) < total_perimeter) then
!             values(1) = 0.0
!             values(2) = 0.0
!             values(3) = total_dPdZ
!         else
!             values(1) = 0.0
!             values(2) = total_dPdZ
!             values(3) = 0.0
!         end if
!     else
!         values(1) = 0.0
!         values(3) = 0.0
!         values(2) = total_dPdZ
!     end if
    
end subroutine

    
!> Compute the gradient of K w.r.t. Z of a cross-section given the depth
subroutine dKdZ(cs, h, value)
    use m_common
    use m_mesh, only: Crosssection
    use m_mesh, only: strickler_type_constant
    use m_mesh, only: strickler_type_powerlaw_h
#ifdef STRICKLER_EINSTEIN
    use m_mesh, only: strickler_type_einstein
#endif
#ifdef AVERAGE_CONVEYANCE
    use m_mesh, only: strickler_type_average_conveyance
#endif
#ifdef DEBORD_FORMULA
    use m_mesh, only: strickler_type_debord
#endif
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Cross-section
    type(Crosssection) :: cs
    !> Water depth
    real(rp), intent(in) :: h
    !> Computed Strickler
    real(rp), intent(out) :: value
    !=LOCAL VARIABLES==================================================================================================!
#ifdef STRICKLER_EINSTEIN
    ! Level
    integer(ip) :: level
    ! Delta h for gradient computation using finite differences
    real(rp), parameter :: dh = 0.1
    ! Computation variables
    real(rp) :: A, B, dA, dB
    ! dPdZ for compound channel
    real(rp), dimension(3) :: dPdZ
    ! Perimeters for compound channel
    real(rp), dimension(3) :: P
    ! Stricklers for gradient computation
    real(rp) :: Kp, Km
#endif

    if (cs%strickler_type_code == strickler_type_constant) then
        value = 0.0_rp
    else if (cs%strickler_type_code == strickler_type_powerlaw_h) then
        value = cs%strickler_params(1) * cs%strickler_params(2) * h**(cs%strickler_params(2) - 1.0_rp)
#ifdef STRICKLER_EINSTEIN
    else if (cs%strickler_type_code == strickler_type_einstein) then
    
        level = cs%level
        call update_level(cs, h+dh)
        call strickler(cs, h+dh, Kp)
        call update_level(cs, h-dh)
        call strickler(cs, h-dh, Km)
        cs%level = level
        value = (Kp - Km) / (2 * dh)
!         if (value > 10.0) value = 10.0
!         if (value < -10.0) value = -10.0
!     
!         ! Compute partial wetted perimeters and derivatives
!         call perimeters_compound_channel(cs, h, P)
!         call dPdZ_compound_channel(cs, h, dPdZ)
!         
!         
!         
!         ! Compute dKdZ
!         A = P(1) + P(2) + P(3)
!         B = P(1) / cs%strickler_params(1)**d3p2 + P(2) / cs%strickler_params(2)**d3p2 + &
!             P(3) / cs%strickler_params(3)**d3p2
!         dA = dPdZ(1) + dPdZ(2) + dPdZ(3)
!         dB = dPdZ(1) / cs%strickler_params(1)**d3p2 + dPdZ(2) / cs%strickler_params(2)**d3p2 + &
!              dPdZ(3) / cs%strickler_params(3)**d3p2
!         value = (dA * B - dB * A) / B**2
#endif
#ifdef AVERAGE_CONVEYANCE
    else if (cs%strickler_type_code == strickler_type_average_conveyance) then
#ifndef CPP_ADJ
        call f90wrap_abort("Subroutine 'dKdZ' must not be called with 'average_conveyance' model")
#endif
#endif
#ifdef DEBORD_FORMULA
    else if (cs%strickler_type_code == strickler_type_debord) then
#ifndef CPP_ADJ
        call f90wrap_abort("Subroutine 'dKdZ' must not be called with 'Debord' model")
#endif
#endif
    end if
    
end subroutine
       

!> Compute the cumulated variable for all cross-sections of the msh
subroutine compute_all_levels_cum(msh)
    use m_common
    use m_mesh, only: Mesh
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Mesh
    type(Mesh), intent(inout) :: msh
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of cross-section
    integer(ip) :: ixs
    
    do ixs = 1, msh%ncs
    
        call compute_levels_cum(msh%cs(ixs))
        
    end do
  
end subroutine

    
!> Update the level of all cross-sections given their respective depths
subroutine update_all_levels(msh, h)
    use m_common
    use m_mesh, only: Mesh
    implicit none
    !=ARGUMENTS========================================================================================================!
    type(Mesh), intent(inout) :: msh
    real(rp), dimension(msh%ncs), intent(in) :: h
    !=LOCAL VARIABLES==================================================================================================!
    ! Index of cross-section
    integer :: ics

    do ics = 1, msh%ncs
        call update_level(msh%cs(ics), h(ics))
    end do
  
end subroutine


! 
!   
! 
! 
! SUBROUTINE SurfaceToHeight( msh,dof )
!    USE m_model
!    USE m_common
!    USE m_mesh
! 
!    implicit none
! 
!    type(Mesh), intent(in)  ::  msh   
!    type( unk ), intent(inout)  ::  dof
! 
!    integer(ip) :: ilevel
!    real(rp) :: square 
!    real(rp) :: s0
!    integer(ip) :: m
!    integer(ip) :: m_max
!    integer(ip) :: findmsurface
!    
!    do ilevel=1,msh%ncs+4
! 
!       !Case 1 : s = 0 (no water) -> h = 0
!       if ((dof%s(ilevel).le.0.0_rp)) then 
!             dof%h(ilevel)=0.0_rp
!       
!       else 
!    
!          ! Compute s0
!          s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
! 
!          !Case 2 : Only rectangular area is wet (s<=s_0 or m =0)
!          if (dof%s(ilevel).le.s0) then
!             dof%h(ilevel)=dof%s(ilevel)/msh%crosssection(ilevel)%surfacewidth(1) ! h= s/w_0
! 
!          !Case 3 : s> s_0 or m>0
!          else
!         
!             m=findMSurface(msh,dof%s(ilevel),ilevel) !return m
!             m_max = size(msh%crosssection(ilevel)%surfacewidth)
! 
!             !Case 3_0 : if m > m_max (rectangular geometry)
!             if (m >= m_max) then
! 
! ! DEPRECATED
! !                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
! !                   (dof%s(ilevel) - msh%crosssection(ilevel)%surfaceareasum(m_max))/msh%crosssection(ilevel)%surfacewidth(m_max)!(s-s_m)/w_m
!                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
!                   (dof%s(ilevel) - s0 - msh%crosssection(ilevel)%trapAreaSum(m_max))/msh%crosssection(ilevel)%surfacewidth(m_max)!(s-s_m)/w_m
!             
!             !Case 3_1 : if w_m+1=w_m (rectangular geometry)
!             else if ((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m)).le.(0.0000001)) then
!                !h=H_m-z_b+hprime_m
!                !hprime_m=(s-s_m)/w_m
! ! DEPRECATED
! !                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
! !                   (dof%s(ilevel)- msh%crosssection(ilevel)%surfaceareasum(m))/msh%crosssection(ilevel)%surfacewidth(m)!(s-s_m)/w_m
!                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
!                   (dof%s(ilevel) - s0 - msh%crosssection(ilevel)%trapAreaSum(m))/msh%crosssection(ilevel)%surfacewidth(m)!(s-s_m)/w_m
! 
!             !Case 3_2 : w_m+1 != w_m (trapezium geometry)
!             else 
!                !h=H_m+hprime_m
!                !hprime_m=-((w_m-square)/(w_{m+1}-w_m))*(H_{m+1}-H_m)
!                ! with square = sqrt( w_m^2+2*((w_{m+1}-w_m)/(H_{m+1}-H_m))*(s-s_m))
! 
! ! DEPRECATED
! !                square=sqrt(msh%crosssection(ilevel)%surfacewidth(m)**2+&                                  !sqrt(w_m^2+
! !                   two*((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m))/ & !2*((w_{m+1}-w_m)/
! !                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m)))* &               !(H_{m+1}-H_m))*
! !                   (dof%s(ilevel)- msh%crosssection(ilevel)%surfaceareasum(m)))                                 !(s-s_m))
!                square=sqrt(msh%crosssection(ilevel)%surfacewidth(m)**2+&                                  !sqrt(w_m^2+
!                   two*((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m))/ & !2*((w_{m+1}-w_m)/
!                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m)))* &               !(H_{m+1}-H_m))*
!                   (dof%s(ilevel) - s0 - msh%crosssection(ilevel)%trapAreaSum(m)))                                 !(s-s_m))
! 
!                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)-&                                !H_m-
!                   (( msh%crosssection(ilevel)%surfacewidth(m)-square)/&                                   !((w_m-square)/
!                   (msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m)))* &   !(w_{m+1}-w_m))*
!                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m))                   !(H_{m+1}-H_m)
!             endif
! 
!          end if
!       end if      
!    end do
! END SUBROUTINE SurfaceToHeight
! 
! !> Update the unknown \f$h_i\f$ from cross section geometry,\f$m_i\f$ and \f$s_i\f$.
! !! \details Algorithm used:\n
! !!  if \f$s_i<=0\f$:
! !!    - \f$h_i=0\f$ 
! !!
! !!  if \f$s_i<=s_i^0\f$:
! !!    - \f$h_i=\frac{s_i}{w_i^0}\f$
! !!
! !!  if \f$s_i>S_i^0\f$:
! !!    - Find \f$m\f$
! !!    - if \f$w_i^{m+1}=w_i^{m}\f$:
! !!       - \f$h_i=H_i^m-z_bi+h'\f$ with \f$h'= \frac{(s_i-s_i^m)}{w_i^m}\f$
! !!    - else: 
! !!       - \f$h_i=H_i^m-z_bi+h'\f$ with \f$h'=-\frac{w_i^m-A}{w_i^{m+1}-w_i^m}(H_i^{m+1}-H_i^m)\f$
! !!       with \f$A = \sqrt{ (w_i^m)^2+2 \frac{w_i^{m+1}-w_i^m}{H_i^{m+1}-H_i^m}(s_i-s_i^m)}\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[inout] dof  Unknowns of the model.
! !! \param[in] crossSection Number of the cross section.
! !! \param[inout] mvector  value of \f$m\f$.
! SUBROUTINE SurfaceToHeightCrossSection( msh,dof,crossSection,mvector)
!    USE m_model
!    USE m_common
!    USE m_mesh
! 
!    implicit none
! 
!    type(Mesh), intent(in)     ::  msh
!    type( unk ), intent(inout)  ::  dof
!    integer(ip), intent(in)     :: crossSection
!    integer(ip), dimension(msh%ncs+4), intent(in)     ::mvector
! 
!    integer(ip) :: ilevel
!    real(rp) :: square
!    real(rp) :: s0
!    integer(ip) :: m
!    integer(ip) :: m_max
!    
!    real(rp) :: StoH
! 
!       ilevel=crossSection
! 
!       !Case 1 : s = 0 (no water) -> h = 0
!       if ((dof%s(ilevel).le.0.0_rp)) then 
!             dof%h(ilevel)=0.0_rp
!       
!       else 
!    
!          ! Compute s0
!          s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
! 
!          !Case 2 : Only rectangular area is wet (s<=s_0 or m =0)
!          if (dof%s(ilevel).le.s0) then
!             dof%h(ilevel)=dof%s(ilevel)/msh%crosssection(ilevel)%surfacewidth(1) ! h= s/w_0
! !             print *, "case2", ilevel, dof%h(ilevel)
! 
!          !Case 3 : s> s_0 or m>0
!          else
!         
!             m=mvector(ilevel) !return m
!             m_max = size(msh%crosssection(ilevel)%surfacewidth)
! 
!             !Case 3_0 : if m>=m_max (rectangular geometry)
!             if (m >= m_max) then
!                
! ! DEPRECATED
! !                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
! !                   (dof%s(ilevel)- msh%crosssection(ilevel)%surfaceareasum(m_max))/msh%crosssection(ilevel)%surfacewidth(m_max)!(s-s_m)/w_m
!                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
!                   (dof%s(ilevel) - s0 - msh%crosssection(ilevel)%trapAreaSum(m_max))/msh%crosssection(ilevel)%surfacewidth(m_max)!(s-s_m)/w_m
!             
! !                print *, "case3_0", ilevel, dof%h(ilevel)
!             !Case 3_1 : if w_m+1=w_m (rectangular geometry)
!             !if (msh%crosssection(ilevel)%surfacewidth(m+1).eq.msh%crosssection(ilevel)%surfacewidth(m)) then
!             else if (abs(msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m)).le.(0.0000001)) then
!                !h=H_m-z_b+hprime_m
!                !hprime_m=(s-s_m)/w_m
!                !dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-msh%crosssection(ilevel)%elevation+&                  !h=H_m-z_b+
! !                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
! !                   (dof%s(ilevel)- msh%crosssection(ilevel)%surfaceareasum(m))/msh%crosssection(ilevel)%surfacewidth(m)!(s-s_m)/w_m
!                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
!                   (dof%s(ilevel) - s0 - msh%crosssection(ilevel)%trapAreaSum(m))/msh%crosssection(ilevel)%surfacewidth(m)!(s-s_m)/w_m
! !                print *, "case3_1", ilevel, dof%h(ilevel)
! 
!             !Case 3_2 : w_m+1 != w_m (trapezium geometry)
!             else 
!                !h=H_m+hprime_m
!                !hprime_m=-((w_m-square)/(w_{m+1}-w_m))*(H_{m+1}-H_m)
!                ! with square = sqrt( w_m^2+2*((w_{m+1}-w_m)/(H_{m+1}-H_m))*(s-s_m))
! 
! ! DEPRECATED
! !                square=sqrt(msh%crosssection(ilevel)%surfacewidth(m)**2+&                                  !sqrt(w_m^2+
! !                   two*((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m))/ & !2*((w_{m+1}-w_m)/
! !                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m)))* &               !(H_{m+1}-H_m))*
! !                   (dof%s(ilevel)- msh%crosssection(ilevel)%surfaceareasum(m)))                                 !(s-s_m))
!                square=sqrt(msh%crosssection(ilevel)%surfacewidth(m)**2+&                                  !sqrt(w_m^2+
!                   two*((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m))/ & !2*((w_{m+1}-w_m)/
!                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m)))* &               !(H_{m+1}-H_m))*
!                   (dof%s(ilevel) - s0 - msh%crosssection(ilevel)%trapAreaSum(m)))                                 !(s-s_m))
! 
!                !dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-msh%crosssection(ilevel)%elevation-&               !H_m-
!                dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)-&               !H_m-
!                   (( msh%crosssection(ilevel)%surfacewidth(m)-square)/&                                   !((w_m-square)/
!                   (msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m)))* &   !(w_{m+1}-w_m))*
!                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m))                   !(H_{m+1}-H_m)
! !                print *, "case3_2", ilevel, dof%h(ilevel), dof%s(ilevel), msh%crosssection(ilevel)%surfacewidth(m), msh%crosssection(ilevel)%surfacewidth(m+1)
! !                print *, "StoH:", StoH(msh, ilevel, dof%S(ilevel))
! !                print *, "bathy:", bathy_cell(ilevel), msh%cs(ilevel)%elevation
! !                print *, "Z:", msh%crosssection(ilevel)%height(:)
! !                print *, "W:", msh%crosssection(ilevel)%surfacewidth(:)
! !                print *, "A:", msh%crosssection(ilevel)%trapAreaSum(:) + s0
! !                print *, "m=", m
!                if (dof%h(ilevel) < -0.1_rp) then
!                   print *, "dof%h(ilevel) < -0.1_rp", dof%h(ilevel), ilevel
!                   stop 1
!                end if
!             endif
! 
!          end if
!       end if
! 
! END SUBROUTINE SurfaceToHeightCrossSection
! 
! !> Determine all value of \f$m\f$ given cross-sections and a wet surface \f$s\f$.
! !! \details Algorithm used:\n
! !!  For each node \f$i\f$:
! !!    - Find \f$m\f$ such that \f$s_i^{m-1}<=s_i<s_i^m\f$ 
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] dof Unknowns of the model.
! !! \param[inout] mvector \f$m\f$ vector to update.
! SUBROUTINE UpdateMVector(msh,dof,mvector)
!    USE m_mesh
!    USE m_model
!    USE m_common
! 
!    implicit none
!    
!    type(Mesh), intent(in)  ::  msh
!    type( unk ), intent(in)  ::  dof
!    integer(ip), dimension(msh%ncs+4),intent(inout)  ::  mvector 
!    
!    real(rp) :: s0
!    integer(ip):: ilevel
! 
!    do ilevel=1,msh%ncs+4
!       mvector(ilevel)=1
! ! DEPRECATED
! !       do while ((msh%crosssection(ilevel)%surfaceareasum(mvector(ilevel))&
! !                .le.dof%s(ilevel)))!.or.(Flood)) ! while s_m < s   
! !          mvector(ilevel)=mvector(ilevel)+1                                               ! m=m+1
! !          if (mvector(ilevel) > size(msh%crosssection(ilevel)%surfaceareasum)) exit
! !       end do
!       s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
!       do while ((msh%crosssection(ilevel)%trapAreaSum(mvector(ilevel)) <= dof%s(ilevel) - s0)) ! while s_m < s   
!          mvector(ilevel)=mvector(ilevel)+1                                               ! m=m+1
!          if (mvector(ilevel) > size(msh%crosssection(ilevel)%trapAreaSum)) exit
!       end do
!       mvector(ilevel)=mvector(ilevel)-1
!    end do
! 
! END SUBROUTINE UpdateMVector
! 
! 
! !> Determine the value of \f$m\f$ given a cross-section (msh + index \f$i\f$ of the cross-section) and a wet surface \f$s\f$.
! !! \details Algorithm used:\n
! !!  For a node \f$i\f$:
! !!    - Find \f$m\f$ such that \f$s_i^{m-1}<=s_i<s_i^m\f$ 
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] surface Wet surface.
! !! \param[in] ilevel Index of the cross section.
! !! \return \f$m\f$ update
! integer(kind=ip) FUNCTION findMSurface(msh,surface,ilevel)
! 
!    USE m_mesh
!    USE m_model
!    implicit none
!    
!    type(Mesh), intent(in)  ::  msh
!    real(rp),intent(in) :: surface
!    integer(ip),intent(in) ::ilevel
!    
!    real(rp) :: s0
!    
!    findMSurface=1
! 
! ! DEPRECATED
! ! TODO : Modifier la condition pour eviter la boucle infinie (verifier que k surfaceareasum(k) existe)
! !    do while ((msh%crosssection(ilevel)%surfaceareasum(findMSurface).le.surface))!.or.(Flood)) ! while s_m < s 
! !         findMSurface=findMSurface+1                                               ! m=m+1    
! !         if (findMSurface > size(msh%crosssection(ilevel)%surfaceareasum)) exit
! !    end do
!    s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
!    do while ((msh%crosssection(ilevel)%trapAreaSum(findMSurface).le.surface - s0)) ! while s_m < s 
!         findMSurface=findMSurface+1                                               ! m=m+1    
!         if (findMSurface > size(msh%crosssection(ilevel)%trapAreaSum)) exit
!    end do
!    findMSurface=findMSurface-1
! END FUNCTION findMSurface
! 
! 
! !> Determine the value of \f$m\f$ given a cross-section (msh + index \f$i\f$ of the cross-section) and a wet surface \f$s\f$.
! !! \details Algorithm used:\n
! !!  For a node \f$i\f$:
! !!    - Find \f$m\f$ such that \f$s_i^{m-1}<=s_i<s_i^m\f$ 
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] dof Unknowns of the model.
! !! \param[in] ilevel Number of the cross section.
! !! \param[inout] mvector \f$m\f$ update
! SUBROUTINE UpdateMVectorElement(msh,dof,ilevel,mvector)
!    USE m_mesh
!    USE m_model
!    USE m_common
! 
!    implicit none
!    
!    type(Mesh), intent(in)  ::  msh
!    type( unk ), intent(in)  ::  dof
!    integer(ip), intent(in)  ::  ilevel
!    integer(ip), dimension(msh%ncs+4),intent(inout)  ::  mvector 
!    
!    real(rp) :: s0
!    
!    mvector(ilevel)=1
! !    do while ((msh%crosssection(ilevel)%surfaceareasum(mvector(ilevel))&
! !       .le.dof%s(ilevel)))!.or.(Flood)) ! while s_m < s   
! !       mvector(ilevel)=mvector(indexcrossSection)+1                                               ! m=m+1
! !       if (mvector(ilevel) > size(msh%crosssection(ilevel)%surfaceareasum)) exit
! !    end do
!    s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
!    do while ((msh%crosssection(ilevel)%trapAreaSum(mvector(ilevel)) <= dof%s(ilevel) - s0)) ! while s_m < s   
!       mvector(ilevel)=mvector(ilevel)+1                                               ! m=m+1
!       if (mvector(ilevel) > size(msh%crosssection(ilevel)%trapAreaSum)) exit
!    end do
!    mvector(ilevel)=mvector(ilevel)-1
! 
! END SUBROUTINE UpdateMVectorElement
! 
! !> Determine the value of \f$m\f$ given a cross-section (msh + index \f$i\f$ of the cross-section) and a water height \f$h\f$.
! !! \details Algorithm used:\n
! !!  For a node \f$i\f$:
! !!    - Find \f$m\f$ such that \f$h_i^{m-1}<=h_i<h_i^m\f$ 
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] height Water height.
! !! \param[in] ilevel Number of the cross section.
! !! \return \f$m\f$ update
! integer(kind=ip) FUNCTION findMHeight(msh,height,ilevel)
!    USE m_mesh
!    USE m_model
! 
!    implicit none
!    
!    type(Mesh), intent(in)  ::  msh
!    real(rp),intent(in) :: height
!    integer(ip),intent(in) ::ilevel
!    
!    real(rp) :: zb
!    
!    findMHeight=1
! 
! ! TODO : Modifier la condition pour eviter la boucle infinie (verifier que k surfaceareasum(k) existe)
!    !zb=msh%crosssection(indexcrossSection)%elevation
!    !zb=bathy_cell(indexcrossSection)
!    do while (msh%crosssection(ilevel)%Height(findMHeight).le.&
!             (height+bathy_cell(ilevel))) ! while H_m < h 
!       findMHeight=findMHeight+1                                                          ! m=m+1
!       if (findMHeight > size(msh%crosssection(ilevel)%Height)) exit
!    end do
!    findMHeight=findMHeight-1
! END FUNCTION findMHeight
! 
! 
! !> Update the unknown \f$h_i\f$ from cross section geometry,\f$m_i\f$ and \f$s_i\f$.
! !! \details Algorithm used:\n
! !!  if \f$s_i<=0\f$:
! !!    - \f$h_i=0\f$ 
! !!
! !!  if \f$s_i<=s_i^0\f$:
! !!    - \f$h_i=\frac{s_i}{w_i^0}\f$
! !!
! !!  if \f$s_i>S_i^0\f$:
! !!    - Find \f$m\f$
! !!    - if \f$w_i^{m+1}=w_i^{m}\f$:
! !!       - \f$h_i=H_i^m-z_bi+h'\f$ with \f$h'= \frac{(s_i-s_i^m)}{w_i^m}\f$
! !!    - else: 
! !!       - \f$h_i=H_i^m-z_bi+h'\f$ with \f$h'=-\frac{w_i^m-A}{w_i^{m+1}-w_i^m}(H_i^{m+1}-H_i^m)\f$
! !!       with \f$A = \sqrt{ (w_i^m)^2+2 \frac{w_i^{m+1}-w_i^m}{H_i^{m+1}-H_i^m}(s_i-s_i^m)}\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] ilevel  Number of the cross section.
! !! \param[in] surface Water surface.
! !! \return Value of the unknow \f$h\f$.
! real (kind=rp) FUNCTION StoH( msh,ilevel, surface )
!    USE m_model
!    USE m_common
!    USE m_mesh
! 
!    implicit none
! 
!    type(Mesh), intent(in)     :: msh
!    integer(ip), intent(in)     :: ilevel
!    real(rp), intent(in)     :: surface
! 
!    real(rp) :: s0
!    real(rp) :: square
!    integer(ip) :: m
!    integer(ip) :: m_max
!    integer(ip) :: findmSurface
!       
! 
!       !Case 1 : s = 0 (no water) -> h = 0
!       if ((surface.le.0.0_rp)) then 
!             StoH=0.0_rp
!       
!       else 
!    
!          ! Compute s0
!          s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
! 
!          !Case 2 : Only rectangular area is wet (s<=s_0 or m =0)
!          if (surface.le.s0) then
!             StoH=surface/msh%crosssection(ilevel)%surfacewidth(1) ! h= s/w_0
! 
!          !Case 3 : s> s_0 or m>0
!          else
!             m=findmSurface(msh,surface,ilevel)!return m
!             m_max = size(msh%crosssection(ilevel)%surfacewidth)
! 
!             !Case 3_0 : if m > m_max (rectangular geometry)
!             if (m >= m_max) then
! 
! ! DEPRECATED
! !                StoH=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
! !                   (surface- msh%crosssection(ilevel)%surfaceareasum(m_max))/msh%crosssection(ilevel)%surfacewidth(m_max)!(s-s_m)/w_m
!                StoH = msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
!                      (surface - s0 - msh%crosssection(ilevel)%trapAreaSum(m_max))/&
!                      msh%crosssection(ilevel)%surfacewidth(m_max)!(s-s_m)/w_m
!             
!             !Case 3_1 : if w_m+1=w_m (rectangular geometry)
!             !if (msh%crosssection(ilevel)%surfacewidth(m+1).eq.msh%crosssection(ilevel)%surfacewidth(m)) then
!             else if (abs(msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m)).le.(0.0000001)) then
!                !h=H_m-z_b+hprime_m
!                !hprime_m=(s-s_m)/w_m
!                !dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-msh%crosssection(ilevel)%elevation+&                  !h=H_m-z_b+
! ! DEPRECATED
! !                StoH=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
! !                   (surface- msh%crosssection(ilevel)%surfaceareasum(m))/msh%crosssection(ilevel)%surfacewidth(m)!(s-s_m)/w_m
!                StoH = msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)+&                                   !h=H_m-z_b+
!                       (surface - s0 - msh%crosssection(ilevel)%trapAreaSum(m))/&
!                       msh%crosssection(ilevel)%surfacewidth(m)!(s-s_m)/w_m
! 
!             !Case 3_2 : w_m+1 != w_m (trapezium geometry)
!             else 
!                !h=H_m+hprime_m
!                !hprime_m=-((w_m-square)/(w_{m+1}-w_m))*(H_{m+1}-H_m)
!                ! with square = sqrt( w_m^2+2*((w_{m+1}-w_m)/(H_{m+1}-H_m))*(s-s_m))
! 
! ! DEPRECATED
! !                square=sqrt(msh%crosssection(ilevel)%surfacewidth(m)**2+&                                  !sqrt(w_m^2+
! !                   two*((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m))/ & !2*((w_{m+1}-w_m)/
! !                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m)))* &               !(H_{m+1}-H_m))*
! !                   (surface- msh%crosssection(ilevel)%surfaceareasum(m)))                                 !(s-s_m))
!                square=sqrt(msh%crosssection(ilevel)%surfacewidth(m)**2+&                                  !sqrt(w_m^2+
!                   two*((msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m))/ & !2*((w_{m+1}-w_m)/
!                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m)))* &               !(H_{m+1}-H_m))*
!                   (surface - s0 - msh%crosssection(ilevel)%trapAreaSum(m)))                                 !(s-s_m))
! 
!                !dof%h(ilevel)=msh%crosssection(ilevel)%height(m)-msh%crosssection(ilevel)%elevation-&               !H_m-
!                StoH=msh%crosssection(ilevel)%height(m)-bathy_cell(ilevel)-&               !H_m-
!                   (( msh%crosssection(ilevel)%surfacewidth(m)-square)/&                                   !((w_m-square)/
!                   (msh%crosssection(ilevel)%surfacewidth(m+1)-msh%crosssection(ilevel)%surfacewidth(m)))* &   !(w_{m+1}-w_m))*
!                   (msh%crosssection(ilevel)%height(m+1)-msh%crosssection(ilevel)%height(m))                   !(H_{m+1}-H_m)
!             endif
! 
!          end if
!       end if
! 
! END FUNCTION StoH
! 
! 
! 
! !> Return the value of the unknown \f$s\f$ from number of the cross section geometry \f$i\f$ and \f$h\f$ value.
! !! \details Algorithm used:\n
! !!  if \f$h<=0\f$:
! !!    - \f$s=0\f$ 
! !!
! !!  if \f$h<=h_i^0\f$:
! !!    - \f$s=h w_i^0\f$
! !!
! !!  if \f$h>h_i^0\f$:
! !!    - Find \f$m\f$
! !!    -\f$s=s_i^m+s'\f$ with \f$s'=\frac{1}{2}(2w_i^m+(\sqrt{(Pe1_i^m)^2-(H_i^{m+1}-H_i^{m})^2}+
! !! \sqrt{(Pe2_i^m)^2-(H_i^{m+1}-H_i^{m})^2})(\frac{h+z_bi-H_i^m}{H_i^{m+1}-H_i^{m}}))(h+z_bi-H_i^m)\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] ilevel  Number of the cross section.
! !! \param[in] height Water height.
! !! \return Value of the unknow \f$s\f$.
! real (kind=rp) FUNCTION HtoS( msh,ilevel, height )
!    USE m_mesh
!    USE m_model
!    implicit none
! 
!    type(Mesh), intent(in)  ::  msh
!    integer(ip),intent(in) :: ilevel
!    real(rp),intent(in) :: height
! 
!    integer(ip) :: findmHeight   
!    integer(ip) :: m
!    integer(ip) :: m_max
!    real(rp) :: s0
! 
!    !Temporary variable
!    real(rp) :: wmp1,wm     !Corresponding to W_{m}
!    real(rp) :: Pe1m,Pe2m   !Corresponding to Pe1_{m} and Pe2_{m}
!    real(rp) :: Hmp1,Hm     !Corresponding to H_{m+1} and H_{m}
!    real(rp) :: zb          !Corresponding to z_b
!    real(rp) :: ainterp     !Interpolation factor
! 
!    HtoS=0.0_rp
! 
!    !Case 1 : height lower than elevation
!    !if (((height-msh%crosssection(ilevel)%elevation).le.0.0_rp)) then
!    if (height.le.0.0_rp) then 
!          HtoS=0.0_rp
!    
!    else 
!       !Case 2 : h< (H_0-z_b)
!       !if (height.le.(msh%crosssection(ilevel)%height(1)-msh%crosssection(ilevel)%elevation)) then
!       if (height.le.(msh%crosssection(ilevel)%height(1)-bathy_cell(ilevel))) then
!       
!          HtoS=(height)*msh%crosssection(ilevel)%surfacewidth(1) !S=h*w_0
! !          print *, "case_2", HtoS
!          
!       else
!       
!          m=findmHeight(msh,height,ilevel)
!          m_max=size(msh%crosssection(ilevel)%height)
!    
!          ! Compute s0
!          s0 = (msh%crosssection(ilevel)%height(1) - bathy_cell(ilevel)) * msh%crosssection(ilevel)%surfacewidth(1)
! !          print *, "S0", msh%crosssection(ilevel)%height(1), bathy_cell(ilevel), msh%crosssection(ilevel)%surfacewidth(1)
! 
!          !Case 3_0 : if m > m_max (rectangular geometry)
!          if (m >= m_max) then
!             Hm =msh%crosssection(ilevel)%height(m_max)
!             wm =msh%crosssection(ilevel)%surfacewidth(m_max)
!             zb =bathy_cell(ilevel)
!             HtoS=msh%crosssection(ilevel)%trapAreaSum(m_max) + s0 + (height+zb-Hm)*wm
! !             print *, "case_3_0", HtoS, msh%crosssection(ilevel)%trapAreaSum(m), s0
!          !Case 3_1 : trapezoidal geometry
!          else
! ! DEPRECATED
! !             wm   =msh%crosssection(ilevel)%surfacewidth(m)
! !             Pe1m =msh%cs(ilevel)%perimeterp1(m+1)
! !             Pe2m =msh%cs(ilevel)%perimeterp2(m+1)
! !             Hm = msh%crosssection(ilevel)%height(m)
! !             Hmp1 = msh%crosssection(ilevel)%height(m+1)
! !             zb = bathy_cell(ilevel)
! ! 
! !             HtoS =  msh%crosssection(ilevel)%trapAreaSum(m) + s0 +   &
! !                  demi*(height+zb-Hm)* &
! !                  ( two*wm + ( sqrt(Pe1m**2-(Hmp1-Hm)**2) + sqrt(Pe2m**2-(Hmp1-Hm)**2) )*  &   
! !                  ((height+zb-Hm)/(Hmp1-Hm)))
! !             print *, "case_3_1", HtoS, msh%crosssection(ilevel)%trapAreaSum(m), s0
!             wm   =msh%crosssection(ilevel)%surfacewidth(m)
!             wmp1 =msh%crosssection(ilevel)%surfacewidth(m+1)
! !             dw = (height+zb-Hm)/(Hmp1-Hm) * (wmp1 + wm)
! !             Pe1m =msh%cs(ilevel)%perimeterp1(m+1)
! !             Pe2m =msh%cs(ilevel)%perimeterp2(m+1)
!             Hm = msh%crosssection(ilevel)%height(m)
!             Hmp1 = msh%crosssection(ilevel)%height(m+1)
!             zb = bathy_cell(ilevel)
! 
!             HtoS =  msh%crosssection(ilevel)%trapAreaSum(m) + s0 +   &
!                     demi*(height+zb-Hm)*((height+zb-Hm)/(Hmp1-Hm)*(wmp1-wm)+two*wm)
! !             print *, "case_3_1", HtoS, msh%crosssection(ilevel)%trapAreaSum(m), s0
! !             print *, "dw=", (height+zb-Hm)/(Hmp1-Hm) * (wmp1 - wm)
! !             print *, "dz=", height+zb-Hm
! !             print *, "w=", 0.5 * (wm + (height+zb-Hm)/(Hmp1-Hm) * (wmp1 - wm) + wm)
! !             print *, "dA=", (height+zb-Hm) * 0.5 * (wm + (height+zb-Hm)/(Hmp1-Hm) * (wmp1 - wm) + wm)
!             
!           end if
!       end if
!    end if
! END FUNCTION HtoS
! 
! 
! 
! !> Return interface value of water surface \f$s\f$ from water elevation \f$H=h+z_b\f$ defined to the interface. 
! !! \details This function computes the value to the equivalent surface to the cell \f$i_l\f$ and \f$i_r\f$
! !! and return the mean value.
! !! \param[in] msh  Mesh of the model.
! !! \param[in] ilevel  Number of the cross section left \f$i_l\f$.
! !! \param[in] icp  Number of the cross section right \f$i_r\f$.
! !! \param[in] height Water elevation  \f$H=h+z_b\f$.
! !! \return Value of the unknow \f$s\f$ to the interface between \f$i_l\f$ and \f$i_r\f$.
! real (kind=rp) FUNCTION HtoSInterface( msh,ilevel,icp, height )
!    USE m_model
!    USE m_common
!    USE m_mesh
!    type(Mesh), intent(in)  ::  msh
!    integer(ip),intent(in) :: ilevel,icp
!    real(rp), intent(in) :: height
! 
!    real(rp) :: HtoS  
!    real(rp) :: slp,srp
! 
!    HtoSInterface=0.0_rp
!    Slp=HtoS(msh,ilevel ,height-bathy_cell(ilevel))
!    Srp=HtoS(msh,icp,height-bathy_cell(icp))
!    HtoSInterface=0.5_rp*(srp+slp)
! END FUNCTION HtoSInterface
! 
! 
! !> Return interface value of water surface \f$s\f$ from water surface defined to the interface. 
! !! \details This function computes the value to the equivalent surface to the cell \f$i_l\f$ and \f$i_r\f$
! !! and return the mean value.
! !! \param[in] msh  Mesh of the model.
! !! \param[in] ilevel  Number of the cross section left \f$i_l\f$.
! !! \param[in] icp  Number of the cross section right \f$i_r\f$.
! !! \param[in] height Water height  \f$h\f$.
! !! \return Value of the unknow \f$s\f$ to the interface between \f$i_l\f$ and \f$i_r\f$.
! real (kind=rp) FUNCTION HtoSInterface_h( msh,ilevel,icp, height )
!    USE m_model
!    USE m_common
!    USE m_mesh
!    type(Mesh), intent(in)  ::  msh
!    integer(ip),intent(in) :: ilevel,icp
!    real(rp), intent(in) :: height
! 
!    real(rp) :: HtoS  
!    real(rp) :: slp,srp
! 
!    Slp=HtoS(msh,ilevel ,height)
!    Srp=HtoS(msh,icp,height)
!    HtoSInterface_h=0.5_rp*(srp+slp)
! 
! END FUNCTION HtoSInterface_h
! 
! 
!> Return the index \f$i\f$ of the cross section the nearest of a point \f$(x,y)\f$.
!! \details This function return the index \f$i\f$ of the cross section the nearest of a point \f$(x,y)\f$.
!! For all points of msh, the function computes the distance with reference points and return the index of the 
!! the point with the lowest distance value.
!! \param[in] msh  Mesh of the model.
!! \param[in] point   \f$(x,y)\f$ point to compare.
!! \return Index \f$i\f$ of the cross section the nearest of the point \f$(x,y)\f$.
FUNCTION search_nearest_cross_section( msh , point ) RESULT( indexi )

    use m_common
    use m_linear_algebra
    USE m_mesh, only: Mesh
    implicit none
    !=ARGUMENTS========================================================================================================!

   type(Mesh), intent(in)  ::  msh
   type(vec2d), intent(in)  ::  point

   integer(ip)  ::  indexi , temp_indexi
   real(rp)    :: xL,yL,xR,yR
   real(rp)    :: distanceTemp,distance

   indexi = -1

   distance=hugem

   do temp_indexi = 3,msh%ncs+2

      xL=point%x
      yL=point%y
      xR=msh%cs(temp_indexi)%coord%x
      yR=msh%cs(temp_indexi)%coord%y
      distanceTemp=sqrt((xR-xL)**2+(yR-yL)**2) !Compute the distance with the point

      
      if (distanceTemp.le.0._rp) then !If point=msh%crosssection(i)%coord 
         indexi=temp_indexi           !   so index =i and exit function
         exit

      else if (distanceTemp.le.distance) then !If distance < distanceTemp
         distance=distanceTemp                ! Get new distanceTemp
         indexi=temp_indexi                   ! Get index i

      endif

   end do

END FUNCTION search_nearest_cross_section
! 
! 
! !> Return the value of \f$w\f$ (water surface width) from msh, water height \f$h\f$ and index of cross-section \f$i\f$.
! !! \details This function return the index \f$i\f$ of the cross section the nearest of a point \f$(x,y)\f$.
! !! For all points of msh, the function computes the distance with reference points and return the index of the 
! !! the point with the lowest distance value.
! !! \details Algorithm used:\n
! !!  if \f$h<=0\f$:
! !!    - \f$w=w_i^0\f$ 
! !!
! !!  if \f$h<=h_i^0\f$:
! !!    - \f$w=w_i^0\f$ 
! !!
! !!  if \f$h>h_i^0\f$:
! !!    - Find \f$m\f$
! !!    - \f$w=w_i^m+\alpha_1(h+z_bi-H_i^m)+\alpha_2(h+z_bi-H_i^m)\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] height Water height \f$h\f$.
! !! \param[in] indexie Indexe \f$i\f$ of the cross section considerated.
! !! \return The value of \f$w\f$ (water surface width)
! Function GetWFromH(msh,height,indexie) result (width)
!    USE m_mesh
!    USE m_model
!    implicit none
! 
!    type(Mesh), intent(in)  ::  msh
!    real(rp), intent (in) :: height
!    integer(ip), intent (in) :: indexie
! 
! 
!    integer(ip) :: findmHeight   
!    integer(ip) :: m
!    integer(ip) :: m_max
!    real(rp) :: width,wm,alpha1,alpha2,Hm,zb
! 
!    !Case 1 : height lower than elevation
!    if (height.le.0.0_rp) then 
!          width=msh%crosssection(indexie)%surfacewidth(1)
!          
!          if (width < 1e-12_rp) then
!             print *, "WRONG WIDTH (case 1)", width
!          end if
!    
!    else 
!       !Case 2 : h< (H_0-z_b)
!       if (height.le.(msh%crosssection(indexie)%height(1)-bathy_cell(indexie))) then
!          width= msh%crosssection(indexie)%surfacewidth(1) !w=w_0
!          
!          if (width < 1e-12_rp) then
!             print *, "WRONG WIDTH (case 2)", width
!          end if
! 
!       !Case 3 : s> s_0 or m>0
!       else
!          m=findmHeight(msh,height,indexie)
!          m_max = size(msh%crosssection(indexie)%surfacewidth)
!          
!          !Case 3_0 : if m > m_max (rectangular geometry)
!          if (m >= m_max) then
!             wm =msh%crosssection(indexie)%surfacewidth(m)
!             width=wm
!             if (width < 1e-12_rp) then
!                 print *, "WRONG WIDTH (case 3-0)", width
!                 print *, "**indexie=", indexie
!             end if
!          !Case 3_1 and 3_2 : if m < m_max
!          else
!             !w=w_{m}+alpha1(H_{m}-h)+ alpha2(H_{m}-h)
!             wm    =msh%crosssection(indexie)%surfacewidth(m)
!             alpha1=msh%cs(indexie)%alpha1(m+1)
!             alpha2=msh%cs(indexie)%alpha2(m+1)
!             Hm    =msh%crosssection(indexie)%height(m)
!             !zb    =msh%cs(indexie)%elevation
!             zb    =bathy_cell(indexie)
!             
!             width=wm+alpha1*(height+zb-Hm)+alpha2*(height+zb-Hm)
!             if (width < 1e-12_rp) then
!                 print *, "WRONG WIDTH (case 3-1)", width
!                 print *, "**indexie=", indexie
!                 print *, "**alpha=", alpha1, alpha2
!                 print *, "**wm,Hm=", wm, Hm
!                 print *, "**height+zb-Hm=", height+zb-Hm
!             end if
!          end if
!       end if
!    end if
! 
! 
! end function GetWFromH
! 
! 
! !> Return the value of perimeter \f$pe\f$ from msh, water height \f$h\f$ and index of cross-section \f$i\f$.
! !! \details 
! !! Algorithm used:\n
! !!  if \f$h<=0\f$:
! !!    - \f$pe=0\f$ 
! !!
! !!  if \f$h<=h_i^0\f$:
! !!    - \f$pe=w_i^0+2h\f$ 
! !!
! !!  if \f$h>h_i^0\f$:
! !!    - Find \f$m\f$
! !!    - \f$pe=pe_i^m+\sqrt{(\frac{w_i^{m+1}}{2}-(\frac{1}{2}w_i^{m}-y_i^{m}))^2+(H_i^{m+1}-H_i^{m})^2}
! !!\frac{(h_i+z_bi-H_i^m)}{H_i^{m+1}-H_i^m}+\sqrt{(\frac{w_i^{m+1}}{2}-(\frac{1}{2}w_i^{m}+y_i^{m}))^2+(H_i^{m+1}-H_i^{m})^2}
! !!\frac{(h_i+z_bi-H_i^m)}{H_i^{m+1}-H_i^m}\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] height Water height \f$h\f$.
! !! \param[in] indexie Indexe \f$i\f$ of the cross section considerated.
! !! \param[in] mvector Vector \f$m\f$.
! !! \return The value of \f$pe\f$ (perimeter)
! Function PerimeterFromH(msh,height,indexie,mvector) result (perimeter)
!    USE m_mesh
!    USE m_model
!    implicit none
! 
!    type(Mesh), intent(in)  ::  msh
!    real(rp), intent (in) :: height
!    integer(ip), intent (in) :: indexie
!    integer(ip), dimension(msh%ncs+4),intent(in)  ::  mvector 
! 
!    real(rp) :: perimeter0
!    integer(ip) :: m, m_max
!    real(rp) :: perimeter,wm,wmp1,ym,Hm,HmP1,zb
! 
!    !Case 1 : height lower than elevation
!    if (height.le.0.0_rp) then 
!          perimeter=0._rp
!    
!    else 
!       !Case 2 : h< (H_0-z_b)
!       if (height.le.(msh%crosssection(indexie)%height(1)-bathy_cell(indexie))) then
! 
!          perimeter= two*(height)+msh%crosssection(indexie)%surfacewidth(1) !pe=w_0+2*h
!          
!       else
!          m=mvector(indexie)
!          m_max=size(msh%crosssection(indexie)%surfacewidth)
!          
!          ! Compute perimeter0
!          perimeter0 = two * (msh%crosssection(indexie)%height(1)-bathy_cell(indexie)) + &
!                       msh%crosssection(indexie)%surfacewidth(1)
!          
!          
!          if (m==m_max) then
!             Hm    =msh%crosssection(indexie)%height(m)
!             zb    =bathy_cell(indexie)
! ! DEPRECATED
! !             perimeter=msh%crosssection(indexie)%perimeterSum(m)
!             perimeter=msh%crosssection(indexie)%trapPerimSum(m) + perimeter0
!             perimeter=perimeter+2.0 * (height+zb-Hm)
!          else
!             wm    =msh%crosssection(indexie)%surfacewidth(m)
!             wmp1  =msh%crosssection(indexie)%surfacewidth(m+1)
!             ym    =msh%crosssection(indexie)%y(m)-msh%crosssection(indexie)%y(m+1)
!             Hm    =msh%crosssection(indexie)%height(m)
!             HmP1  =msh%crosssection(indexie)%height(m+1)
!             !zb    =msh%crosssection(indexie)%elevation
!             zb    =bathy_cell(indexie)
! ! DEPRECATED
!             !             perimeter=msh%crosssection(indexie)%perimeterSum(m)
! 
!             perimeter=msh%crosssection(indexie)%trapPerimSum(m) + perimeter0
! 
!             !See doc Pe(h)
!             perimeter=perimeter+&
!                       sqrt((demi*wmp1-(demi*wm-ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))+&
!                       sqrt((demi*wmp1-(demi*wm+ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))
!          end if
!          !<NOADJ
! !           if (perimeter < 1e-12) then
! ! 
! !             print *, "perimeter <= 0"
! !             print *, indexie, m, msh%crosssection(indexie)%trapPerimSum(m), perimeter0
! !             print *, height + zb, Hm
! !             print *, msh%crosssection(indexie)%height(1:m)
! !             read(*,*)
! !             
! !           end if
!          !>NOADJ
!          
!       end if
!    end if
! 
! end function PerimeterFromH
! 
! 
! !> Return the value of perimeter \f$pe\f$ from msh, water height \f$h\f$ and index of cross-section \f$i\f$.
! !! \details 
! !! Algorithm used:\n
! !!  if \f$h<=0\f$:
! !!    - \f$pe=0\f$ 
! !!
! !!  if \f$h<=h_i^0\f$:
! !!    - \f$pe=w_i^0+2h\f$ 
! !!
! !!  if \f$h>h_i^0\f$:
! !!    - Find \f$m\f$
! !!    - \f$pe=pe_i^m+\sqrt{(\frac{w_i^{m+1}}{2}-(\frac{1}{2}w_i^{m}-y_i^{m}))^2+(H_i^{m+1}-H_i^{m})^2}
! !!\frac{(h_i+z_bi-H_i^m)}{H_i^{m+1}-H_i^m}+\sqrt{(\frac{w_i^{m+1}}{2}-(\frac{1}{2}w_i^{m}+y_i^{m}))^2+(H_i^{m+1}-H_i^{m})^2}
! !!\frac{(h_i+z_bi-H_i^m)}{H_i^{m+1}-H_i^m}\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] height Water height \f$h\f$.
! !! \param[in] indexie Indexe \f$i\f$ of the cross section considerated.
! !! \param[in] mvector Vector \f$m\f$.
! !! \return The value of \f$pe\f$ (perimeter)
! subroutine Perimeters3FromH(msh,height,indexie,mvector, perimeters)
!    USE m_mesh
!    USE m_model
!    implicit none
! 
!    type(Mesh), intent(in)  ::  msh
!    real(rp), intent (in) :: height
!    integer(ip), intent (in) :: indexie
!    integer(ip), dimension(msh%ncs+4),intent(in)  ::  mvector 
! 
!    real(rp) :: perimeter0
!    integer(ip) :: m, m_max
!    real(rp) :: perimeter,wm,wmp1,ym,Hm,HmP1,zb
!    real(rp), dimension(3) :: perimeters
! 
!    !Case 1 : height lower than elevation
!    if (height.le.0.0_rp) then 
!          perimeters(:)=0._rp
!    
!    else 
!       !Case 2 : h< (H_0-z_b)
!       if (height.le.(msh%crosssection(indexie)%height(1)-bathy_cell(indexie))) then
! 
!          perimeters(1)= two*(height)+msh%crosssection(indexie)%surfacewidth(1) !pe=w_0+2*h
!          perimeters(2:3) = 0.0
!          
!       else
!          m=mvector(indexie)
!          m_max=size(msh%crosssection(indexie)%surfacewidth)
!          
!          ! Compute perimeter0
!          perimeter0 = two * (msh%crosssection(indexie)%height(1)-bathy_cell(indexie)) + &
!                       msh%crosssection(indexie)%surfacewidth(1)
!          
!          
!          if (m==m_max) then
!             Hm    =msh%crosssection(indexie)%height(m)
!             zb    =bathy_cell(indexie)
! ! DEPRECATED
! !             perimeter=msh%crosssection(indexie)%perimeterSum(m)
! !             perimeter=msh%crosssection(indexie)%trapPerimSum(m) + perimeter0
! !             perimeter=perimeter+2.0 * (height+zb-Hm)
!             perimeters(1) = perimeter0 + msh%crosssection(indexie)%trapPerimSum(2)
!             perimeters(2) = msh%crosssection(indexie)%trapPerimSum(3) - msh%crosssection(indexie)%trapPerimSum(2)
!             perimeters(3) = msh%crosssection(indexie)%trapPerimSum(4) - msh%crosssection(indexie)%trapPerimSum(3) + &
!                             2.0 * (height+zb-Hm)
!             
!          else
!          
!             wm    =msh%crosssection(indexie)%surfacewidth(m)
!             wmp1  =msh%crosssection(indexie)%surfacewidth(m+1)
!             ym    =msh%crosssection(indexie)%y(m)-msh%crosssection(indexie)%y(m+1)
!             Hm    =msh%crosssection(indexie)%height(m)
!             HmP1  =msh%crosssection(indexie)%height(m+1)
!             zb    =bathy_cell(indexie)
!             
!             if (m == 1) then
!                perimeters(1) = perimeter0 + msh%crosssection(indexie)%trapPerimSum(1) + &
!                                sqrt((demi*wmp1-(demi*wm-ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))+&
!                                sqrt((demi*wmp1-(demi*wm+ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))
!             else if (m == 2) then
!             
!                perimeters(1) = perimeter0 + msh%crosssection(indexie)%trapPerimSum(2)
!                perimeters(2) = sqrt((demi*wmp1-(demi*wm-ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))+&
!                                sqrt((demi*wmp1-(demi*wm+ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))
!                perimeters(3) = 0.0_rp
!                
!             else
!             
!                perimeters(1) = perimeter0 + msh%crosssection(indexie)%trapPerimSum(2)
!                perimeters(2) = msh%crosssection(indexie)%trapPerimSum(3) - msh%crosssection(indexie)%trapPerimSum(2)
!                perimeters(3) = sqrt((demi*wmp1-(demi*wm-ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))+&
!                                sqrt((demi*wmp1-(demi*wm+ym))**2+(Hmp1-hm)**2)*((height+zb-Hm)/(Hmp1-Hm))
!                                
!             end if
!             
!          end if
!          
!       end if
!    end if
! 
! end subroutine Perimeters3FromH
! 
! 
! !> Return the value of \f$\frac{\partial w}{\partial Z}\f$ from msh, water height \f$h\f$ and index of cross-section \f$i\f$.
! !! \details 
! !! \details Algorithm used: TO DO\n
! !!  if \f$h<=0\f$:
! !!    - \f$w=w_i^0\f$ 
! !!
! !!  if \f$h<=h_i^0\f$:
! !!    - \f$w=w_i^0\f$ 
! !!
! !!  if \f$h>h_i^0\f$:
! !!    - Find \f$m\f$
! !!    - \f$w=w_i^m+\alpha_1(h+z_bi-H_i^m)+\alpha_2(h+z_bi-H_i^m)\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] height Water height \f$h\f$.
! !! \param[in] indexie Indexe \f$i\f$ of the cross section considerated.
! !! \return The value of \f$\frac{\partial w}{\partial Z}\f$ (derivative of the width with respect to the elevation)
! Function dWdZFromH(msh,height,indexie) result (dWdZ)
!   USE m_mesh
!   USE m_model
!   implicit none
! 
!   type(Mesh), intent(in)  ::  msh
!   real(rp), intent (in) :: height
!   integer(ip), intent (in) :: indexie
!   
!   integer(ip) :: findmHeight   
!   integer(ip) :: m
!   integer(ip) :: m_max
!   real(rp) :: alpha1,alpha2
!   real(rp) :: dWdZ
! 
!   !Case 1 : h <= H0
!   if (height.le.(msh%crosssection(indexie)%height(1)-bathy_cell(indexie))) then 
!      dWdZ = 0._rp
!      
!   !Case 2 : h > H0
!   else 
!      m=findmHeight(msh,height,indexie)
!      m_max = size(msh%crosssection(indexie)%surfacewidth)
!      
!      !Case 2_0 : if m > m_max (rectangular geometry)
!      if (m >= m_max) then
!         dWdZ = 0._rp
!      !Case 2_1
!      else
!         alpha1=msh%cs(indexie)%alpha1(m+1)
!         alpha2=msh%cs(indexie)%alpha2(m+1)
!         dWdZ = alpha1 + alpha2
!      end if
! 
!   end if
! 
!  end function dWdZFromH
! 
! 
! 
! 
! !> Return the value of \f$\frac{\partial pe}{\partial Z}\f$ from msh, water height \f$h\f$ and index of cross-section \f$i\f$.
! !! \details 
! !! Algorithm used:\n TODO
! !!  if \f$h<=0\f$:
! !!    - \f$pe=0\f$ 
! !!
! !!  if \f$h<=h_i^0\f$:
! !!    - \f$pe=w_i^0+2h\f$ 
! !!
! !!  if \f$h>h_i^0\f$:
! !!    - Find \f$m\f$
! !!    - \f$pe=pe_i^m+\sqrt{(\frac{w_i^{m+1}}{2}-(\frac{1}{2}w_i^{m}-y_i^{m}))^2+(H_i^{m+1}-H_i^{m})^2}
! !!\frac{(h_i+z_bi-H_i^m)}{H_i^{m+1}-H_i^m}+\sqrt{(\frac{w_i^{m+1}}{2}-(\frac{1}{2}w_i^{m}+y_i^{m}))^2+(H_i^{m+1}-H_i^{m})^2}
! !!\frac{(h_i+z_bi-H_i^m)}{H_i^{m+1}-H_i^m}\f$
! !!
! !! \param[in] msh  Mesh of the model.
! !! \param[in] height Water height \f$h\f$.
! !! \param[in] indexie Indexe \f$i\f$ of the cross section considerated.
! !! \param[in] mvector Vector \f$m\f$.
! !! \return The value of \f$pe\f$ (perimeter)
! Function dPdZFromH(msh,height,indexie,mvector) result (dPdZ)
!    USE m_mesh
!    USE m_model
!    implicit none
! 
!    type(Mesh), intent(in)  ::  msh
!    real(rp), intent (in) :: height
!    integer(ip), intent (in) :: indexie
!    integer(ip), dimension(msh%ncs+4),intent(in)  ::  mvector 
! 
! 
!    integer(ip) :: m, m_max
!    real(rp) :: dPdZ,wm,wmp1,ym,Hm,HmP1
! 
!    !Case 1 : height lower than elevation
!    if (height.le.0.0_rp) then 
!          dPdZ = 0._rp
!    
!    else 
!       !Case 2 : h< (H_0-z_b)
!       if (height.lt.(msh%crosssection(indexie)%height(1)-bathy_cell(indexie))) then
!          dPdZ = 2.0_rp
!          
!       else
!          m=mvector(indexie)
!          m_max=size(msh%crosssection(indexie)%surfacewidth)
!          
!          if (m==m_max) then
!             dPdZ = 2.0_rp
!          else
!             wm    =msh%crosssection(indexie)%surfacewidth(m)
!             wmp1  =msh%crosssection(indexie)%surfacewidth(m+1)
!             ym    =msh%crosssection(indexie)%y(m)-msh%crosssection(indexie)%y(m+1)
!             Hm    =msh%crosssection(indexie)%height(m)
!             HmP1  =msh%crosssection(indexie)%height(m+1)
! 
!             !See doc Pe(h)
!             dPdZ = two * sqrt((demi*wmp1-(demi*wm-ym))**2+(Hmp1-hm)**2) / (Hmp1-Hm)
!                       
!          end if
!       end if
!    end if
! 
! 
!  end function dPdZFromH
! 
! 
! 
!  
! 
! !> Return the distance between two points \f$(x_l,y_l)\f$ and \f$(x_r,y_r)\f$.
! !! \details 
! !! Compute the distance \f$d\f$  between this two points ( \f$d=\sqrt{(x_r-x_l)^2+(y_r-y_l)^2}\f$)
! !! \param[in] xl Coordinate \f$x\f$ of the first point.
! !! \param[in] yl Coordinate \f$y\f$ of the first point.
! !! \param[in] xr Coordinate \f$x\f$ of the second point.
! !! \param[in] yr Coordinate \f$y\f$ of the second point.
! !! \return Distance between the two points \f$(x_l,y_l)\f$ and \f$(x_r,y_r)\f$.
! Function distanceBetweenToPoint(xL,yL,xR,yR) result (distance)
!    USE m_mesh
!    implicit none
! 
!    real(rp), intent (in) :: xL,yL,xR,yR
!    real(rp) :: distance 
!    distance=sqrt((xR-xL)**2+(yR-yL)**2)
! end function distanceBetweenToPoint
! 
! 
! 
