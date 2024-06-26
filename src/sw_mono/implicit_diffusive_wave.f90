!======================================================================================================================!
!
!                    DassFlow1D Version 2.1
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA - CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software aiming at simulating geophysical free surface flows.
!  It is designed for Variational Data Assimilation (4D-var) and sensitivity analyses. 
! Inverse capabilities are based on the adjoint code which is generated by 
! a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures:
!    - DassFlow 2DShallow (shallow water equations in (h,q), finite volumes), i.e. the present code.
!    - DassFlow3D  (non-linear Stokes equations, finite elements, mobile geometries, ALE).
!    - DassFlow 1D (shallow water equations in (S,Q), finite volumes or finite differences), river hydraulics.
!  Please consult the DassFlow webpage for more details: http://www.math.univ-toulouse.fr/DassFlow
!
! You have used DassFlow in an article, a scientific document etc ?  How to cite us ? Please consult the webpage.
! 
!  Many people have contributed to the DassFlow developments from the initial version to the latest ones.
!  Current main developers or scientific contributers are:
!               P. Brisset (CNES & Mathematics Institute of Toulouse & INSA Toulouse)
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT)
!               P.-A. Garambois (INSA & ICUBE Strasbourg)
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!  and former other developers (R. Madec, M. Honnorat and J. Marin).
!
!  Scientific contact : jerome.monnier@insa-toulouse.fr
!  Technical  contacts : frederic.couderc@math.univ-toulouse.fr, pierre.brisset@insa-toulouse.fr
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
!> \file implicit_diffusive_wave.f90
!! \brief This file includes the computation with implicit diffusive wave equation methode without MUSCL and slope limiter
!! \details For more details, see developper documentation. This file contains only implicit_diffusive_wave subroutine


!> Subroutine of the computation implicit diffusive wave equation method without MUSCL and without slope limiter. 
!!
!! \details More details given in developper documentation
!! \param[in]  dof Unknowns of the model.
!! \param[in]    msh mshof the model.

SUBROUTINE implicit_diffusive_wave(mdl, msh, large_grid, imp, dof, status)

   USE m_common
   USE m_mesh
   !USE m_time_screen          !NOADJ
   use m_sw_mono, only: ImplicitMatrix
   use m_sw_mono, only: Model
   use m_sw_mono, only: Unknowns
   USE m_linear_solver        
   !USE m_user_data
   USE m_numeric
   implicit none
   
!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE(Mesh), intent(inout)  ::  msh
   TYPE(Mesh), intent(inout)  ::  large_grid
   TYPE(Model), intent(in) :: mdl
   TYPE(ImplicitMatrix), intent(inout) :: imp
   TYPE(Unknowns), intent(inout) :: dof
   integer(ip) :: status

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!
   integer(ip) :: info

   ! Element of the tridiagonal matrice to be invert in the Newton-Raphson equation
   real(rp), dimension(msh%ncs) :: DAL  ! Lower diagonal
   real(rp), dimension(msh%ncs) :: DAU  ! Upper diagonal
   real(rp), dimension(msh%ncs) :: DA   ! Diagonal

   ! Right hand side of the Newton Raphson equation
   real(rp), dimension(msh%ncs) :: R  ! Right hand side of the equation
   real(rp) :: B, C                       ! Free surface height and bathymetry terms
   
   ! Physical variables 
   real(rp), dimension(msh%ncs)  :: v  ! linearized velocity
   real(rp), dimension(msh%ncs)  :: Hk ! linearized free surface height 
   real(rp), dimension(large_grid%ncs)  :: Hk_lg ! linearized free surface height on the large grid 

   ! Geometry 
   real(rp) :: w0u, whu, wbu               ! resp. bed, free surface and average upstream width 
   real(rp) :: w0, wh, wb                  ! resp. bed, free surface and average centered width
   real(rp) :: w0d, whd, wbd               ! resp. bed, free surface and average downstream width
   real(rp) :: whi, wbi                    ! resp. free surface and average initial width 
   real(rp) :: dtauwb, dywb                ! time and space derivative of the average width
   real(rp) :: Rh                          ! hydraulic radius
   real(rp) :: pw                          ! wet perimeter
   
   ! Mesh
   real(rp) :: dxie,dxiem,dxb           ! mean space step of the mesh
   real(rp) :: dyb                      ! mean space step of the large grid
   real(rp) :: dtau                     ! time step
   integer(rp) :: ibc_us,ibc_ds         ! upstream and downstream boundary conditions indexes
   integer(rp) :: iseg                  ! segment index

   ! Physical parameters
   real(rp) :: dxHk,dxxHk   ! slope of the free surface and second order derivative
   real(rp) :: gam        ! puissance de la relation de Manning-Strickler
   real(rp) :: aa,bb,cc,phi,psi ! coefficient of the equation
   real(rp) :: K ! coefficient de la loi de Manning-Strickler
   integer(rp) :: ie ! node index 
   integer(rp) :: ie_large_grid ! node index
   integer(rp) :: ie_us,ie_ds ! index of the up and downstream node
   integer(rp) :: ratio_geom

   ! Newton-Raphson algorithm criteria
   real(rp) :: norm_err    ! norme of the error for the convergence criteria 
   real(rp) :: tol         ! tolerance of the convergence criteria
   integer(ip) :: it_sol   ! Newton's method current iteration   


   !======================================================================================================================!
   !  Begin Subroutine
   !======================================================================================================================!


   !===================================================================================================================!
   !  Newton s method for nonlinear systems
   !===================================================================================================================!
   
   ! initialization of Newton-Raphson algorithm criteria
   norm_err = 1.0_rp
   tol = 1.0e-9_rp
   it_sol = 0
   ie_large_grid = 0
   gam = 5.0/3.0
   if (mdl%scale_idw == "large") then
      ratio_geom = 100
   else 
      ratio_geom = 1
   end if 
   dtau = mdl%dt / ratio_geom

   ! Initialization of the boundary condition
   iseg=1 ! A MODIFIER POUR LE RESEAU
   ie_us = msh%seg(iseg)%first_cs
   ie_ds = msh%seg(iseg)%last_cs

   ! upstream boundary condition
   if (msh%seg(iseg)%us_bc > 0) then
      ibc_us = msh%seg(iseg)%us_bc
      ! Dirichlet boundary condition on the elevation
      if (mdl%bc(ibc_us)%id == "elevation") then 
	 !write(6,*)"test upstream boundary condition",linear_interp(mdl%bc(ibc_us)%ts%t, mdl%bc(ibc_us)%ts%y, mdl%tc) - msh%cs(ie_us)%bathy
         dof%h(ie_us) = linear_interp(mdl%bc(ibc_us)%ts%t, mdl%bc(ibc_us)%ts%y, mdl%tc) - msh%cs(ie_us)%bathy
      ! Dirichlet boundary condition on the discharge
      else if (mdl%bc(ibc_us)%id == "discharge") then
         dof%Q(ie_us) = linear_interp(mdl%bc(ibc_us)%ts%t, mdl%bc(ibc_us)%ts%y, mdl%tc)
      end if 
   ! RAJOUTER UN ELSE POUR LA RESOLUTION SUR LE RESEAU
   end if 

   ! downstream boundary condition
   if (msh%seg(iseg)%ds_bc > 0) then
      ibc_ds = msh%seg(iseg)%ds_bc
      ! Dirichlet boundary condition on the elevation
      if (mdl%bc(ibc_ds)%id == "elevation") then
	 !write(6,*)"test downstream boundary condition",linear_interp(mdl%bc(ibc_ds)%ts%t, mdl%bc(ibc_ds)%ts%y, mdl%tc) - msh%cs(ie_ds)%bathy
         dof%h(ie_ds) = linear_interp(mdl%bc(ibc_ds)%ts%t, mdl%bc(ibc_ds)%ts%y, mdl%tc) - msh%cs(ie_ds)%bathy
      ! Dirichlet boundary condition on the discharge
      else if (mdl%bc(ibc_ds)%id == "discharge") then
         dof%Q(ie_ds) = linear_interp(mdl%bc(ibc_ds)%ts%t, mdl%bc(ibc_ds)%ts%y, mdl%tc)
      end if
   ! RAJOUTER UN ELSE POUR LA RESOLUTION SUR LE RESEAU
   end if

   ! setting the first guess and the space step 
   Hk(ie_us-2:ie_ds+2) = dof%h(ie_us-2:ie_ds+2) + msh%cs(ie_us-2:ie_ds+2)%bathy
         
   ! loop for the Newton Raphson algorithm
   do while ( it_sol <= 100 .and. norm_err > tol )

      if (mdl%scale_idw == "large") then
         ie_large_grid = large_grid%seg(iseg)%first_cs
         do ie = ie_us,ie_ds 
            if ( abs(msh%cs(ie)%x - large_grid%cs(ie_large_grid)%x) < 1.0e-6_rp ) then 
               Hk_lg(ie_large_grid) = Hk(ie)
               ie_large_grid = ie_large_grid + 1
            end if 
          end do
          ie_large_grid = large_grid%seg(iseg)%first_cs
      end if
      

      ! loop inside the domain
      do ie = ie_us+1,ie_ds-1

         ! mean space step 
         dxie = msh%cs(ie+1)%deltademi / ratio_geom
         dxiem = msh%cs(ie)%deltademi / ratio_geom
         !write(6,*) "test space step", dxie, dxiem
         dxb = ( dxie + dxiem ) *0.5 !( dy(ie) + dy(ie-1) ) *0.5

         ! slope and second order derivative of the free surface height
         dxHk = ( Hk(ie+1) - Hk(ie-1) ) / (2*dxb)
         if (abs(dxHk) <= 1.0e-6_rp) then 
            write(6,*) "the slope is freaking low",dxHk
            dxHk = - 1.0e-6_rp 
         end if 
         dxxHk = Hk(ie+1) / ( dxie*dxb )  &
               - ( 2 * Hk(ie) ) / ( dxie*dxiem ) &
               + Hk(ie-1) / ( dxiem*dxb ) 

         ! width variations
	 if (mdl%scale_idw == "large") then
            
            ! update the reach index 
            if (abs(msh%cs(ie)%x-large_grid%cs(ie_large_grid)%x)>abs(msh%cs(ie)%x-large_grid%cs(ie_large_grid+1)%x)) then 
               ie_large_grid = ie_large_grid + 1 
            end if 
            
            ! upstream widths
            call width(large_grid%cs(ie_large_grid+1), Hk_lg(ie_large_grid+1) - large_grid%cs(ie_large_grid+1)%bathy, whu)
            !call width(large_grid%cs(ie_large_grid+1), zero, w0u)
            w0u = large_grid%cs(ie_large_grid+1)%level_widths(1)
            wbu = ( w0u + w0u ) * 0.5

            ! centered widths
            call width(large_grid%cs(ie_large_grid), Hk_lg(ie_large_grid) - large_grid%cs(ie_large_grid)%bathy, wh)
            !call width(large_grid%cs(ie_large_grid), zero, w0)
            w0 = large_grid%cs(ie_large_grid)%level_widths(1)
            wb = ( w0 + w0 ) * 0.5

            ! downstream widths
            call width(large_grid%cs(ie_large_grid-1), Hk_lg(ie_large_grid-1) - large_grid%cs(ie_large_grid-1)%bathy, whd)
            !call width(large_grid%cs(ie_large_grid-1), zero, w0d)
            w0d = large_grid%cs(ie_large_grid-1)%level_widths(1)
            wbd = ( w0d + w0d ) * 0.5

            ! initial widths
            call width(msh%cs(ie), dof%h(ie), whi)
            wbi = ( w0 + w0 ) * 0.5

            dyb = ( large_grid%cs(ie_large_grid+1)%deltademi + large_grid%cs(ie_large_grid)%deltademi ) *0.5 / ratio_geom
            dtauwb = ( wb - wbi ) / dtau
            dywb = ( wbu - wbd ) / (2*dyb)

         else if (mdl%scale_idw == "low") then
            wbu = 0.0_rp
            whu = 0.0_rp
            w0u = 0.0_rp

            wb = 1.0_rp
            wh = 0.0_rp
            w0 = 0.0_rp

            wbd = 0.0_rp
            whd = 0.0_rp
            w0d = 0.0_rp

            wbi = 0.0_rp
            whi = 0.0_rp

            dtauwb = 0.0_rp
            dywb = 0.0_rp
	 end if 
         !write(6,*) "test largeur",wb,dtauwb,dywb
         ! compute the velocity
         call perimeter(msh%cs(ie), Hk(ie) - msh%cs(ie)%bathy, pw)
         Rh = htoA(msh%cs(ie), Hk(ie) - msh%cs(ie)%bathy) / pw
         call strickler(msh%cs(ie), Hk(ie) - msh%cs(ie)%bathy, K)
         K = K/(ratio_geom**0.5)
         !write(6,*) "test strickler", K
         v(ie)= K * abs(dxHk)**0.5 * Rh**(gam-1) 
         !if (v(ie) >= 1.0_rp) then
         !    write(6,*) "mayday"
         !    write(6,*) "test velocity 0",K, abs(dxHk)**0.5, Rh**(gam-1)
         !    write(6,*) "test velocity 0",K, v(ie)
         !end if 
         ! the coefficient of the equation 
         aa = 0.5 * v(ie) * ( Hk(ie) - msh%cs(ie)%bathy ) / abs(dxHk)
         bb = gam * v(ie)
         phi = gam * v(ie) - ( 0.5 * v(ie) * ( Hk(ie) - msh%cs(ie)%bathy ) * dxxHk / (dxHk **2 ) )
          
         cc = ( dtauwb + gam*v(ie)*dywb ) / wb
         psi = (1/dtau) + ( ( dtauwb + gam*v(ie)*dywb ) / wb ) - ( 0.5 * v(ie) * dxxHk / dxHk )
         !write(6,*) "test output", K , abs(dxHk)**0.5 , Rh**(gam-1) 

         ! the matrice A on the domain 
         DAL(ie) = - ( aa/(dxiem*dxb) ) - ( phi / (2*dxb) ) ! the lower diagonal  
         DA(ie) = 2*aa / (dxiem*dxie ) + psi              ! the diagonal 
         DAU(ie) = - ( aa/(dxie*dxb) ) + ( phi / (2*dxb) )   ! the upper diagonal 

         ! product of the matrice B with the vector H on the domain 
         B = ( ( aa/(dxie*dxb) ) - ( bb / (2*dxb) ) ) * ( Hk(ie+1) )   &    ! upperdiagonal
           + ( - 2*aa / (dxiem*dxie ) - cc -( 1 / dtau) ) * ( Hk(ie) )   & ! diagonal
           + ( ( aa/(dxiem*dxb) ) + ( bb / (2*dxb) ) ) * ( Hk(ie-1) )      ! lowerdiagonal

         ! product of the matrice C with the vector Z on the domain
         C =  ( 0.5 * bb / dxb ) * msh%cs(ie+1)%bathy &  ! upperdiagonal
             + cc * msh%cs(ie)%bathy                   & ! diagonal
             - ( 0.5 * bb / dxb ) * msh%cs(ie-1)%bathy   ! lowerdiagonal

         ! the right hand side
         R(ie) = (1/dtau)*(dof%h(ie)+msh%cs(ie)%bathy) + B + C

      end do
      !pause

      ! Boundary conditions

      ! upstream boundary 
      ! Dirichlet boundary condition on the elevation
      if (mdl%bc(ibc_us)%id == "elevation") then
         !write(6,*) "test bord us ",Hk(ie_us),bc_in - Hk(ie_us)
         DA(ie_us) = 1
         DAU(ie_us) = 0
         R(ie_us) = 0
      ! Dirichlet boundary condition on the discharge
      else if (mdl%bc(ibc_us)%id == "discharge") then 
         call perimeter(msh%cs(ie_us), Hk(ie_us) - msh%cs(ie_us)%bathy, pw)
         call strickler(msh%cs(ie_us), Hk(ie_us) - msh%cs(ie_us)%bathy, K)
         Rh = htoA(msh%cs(ie_us), Hk(ie_us) - msh%cs(ie_us)%bathy) / pw
         Hk( ie_us - 1 ) = Hk( ie_us ) + Hk( ie_us ) - Hk( ie_us + 1 )
         !Hk(ie_us-1) = Hk( ie_us ) + msh%cs( ie_us )%bathy - msh%cs( ie_us + 1 )%bathy
         dxie = msh%cs(ie_us+1)%deltademi 
         dxHk = ( Hk(ie_us+1) - Hk(ie_us-1) ) / ( 2*dxie)
         DA(ie_us) = ( 1 / ( Hk(ie_us) - msh%cs( ie_us )%bathy ) ) + ( 1 / ( 2 * abs(dxHk) * dxie) ) 
         DAU(ie_us) = - 1 / (2 * abs(dxHk) * dxie)
         R(ie_us) = ( dof%Q(ie_us) / ( K * abs(dxHk)**(0.5) * pw * Rh**(gam) ) ) - 1 
         !write(6,*) "test bord us ", R(ie_us)
      ! Neuman boundary condition
      else if (mdl%bc(ibc_us)%id == "normal_depth") then 
         DA(ie_us) = 1
         DAU(ie_us) = -1
         R(ie_us) = 0
      end if
      
      ! downstream boundary
      ! Dirichlet boundary condition on the elevation
      if (mdl%bc(ibc_ds)%id == "elevation") then
         DA(ie_ds) = 1
         DAL(ie_ds) = 0  
         R(ie_ds) = 0 
      ! Dirichlet boundary condition on the discharge
      else if (mdl%bc(ibc_ds)%id == "discharge") then 
         call perimeter(msh%cs(ie_ds), Hk(ie_ds) - msh%cs(ie_ds)%bathy, pw)
         call strickler(msh%cs(ie_ds), Hk(ie_ds) - msh%cs(ie_ds)%bathy, K)
         Rh = htoA(msh%cs(ie_ds), Hk(ie_ds) - msh%cs(ie_ds)%bathy) / pw
         Hk( ie_ds + 1 ) = Hk( ie_ds ) + Hk( ie_ds ) - Hk( ie_ds - 1 )
         dxie = msh%cs(ie_ds)%deltademi
         dxHk = ( Hk(ie_ds+1) - Hk(ie_ds-1) ) / ( 2 * dxie )
         DA(ie_ds) = ( 1 / ( Hk(ie_ds) - msh%cs( ie_ds )%bathy ) ) + ( 1 / ( 2 * abs(dxHk) * dxie ) ) 
         DAU(ie_ds) = - 1 / (2 * abs(dxHk) * dxie)
         R(ie_ds) = dof%Q(ie_ds) / ( K * abs(dxHk)**(0.5) * pw * Rh**(gam) ) - 1  
      ! Neuman boundary condition
      else if (mdl%bc(ibc_ds)%id == "normal_depth") then
         Hk( ie_ds + 1 ) = Hk( ie_ds ) + Hk( ie_ds ) - Hk( ie_ds - 1 )
         dxie = msh%cs(ie_ds)%deltademi / ratio_geom
         dxHk = ( Hk(ie_ds+1) - Hk(ie_ds-1) ) / ( 2 * dxie)
         DA(ie_ds) = 1         
         DAL(ie_ds) = -1
         R(ie_ds) = msh%cs( ie_ds )%bathy-msh%cs( ie_ds -1)%bathy-(dxHk * dxie) 
         !write(6,*) "test bord ds ", R(ie_ds)
      end if 
      

      ! update the of the stopping criteria, the maximum of the residu
      norm_err  =  sqrt(sum( R(ie_us:ie_ds)**2))

      ! Call of linear solver (Solve AX=R with A tridiagonal matrix)
      call DGTSV( msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs + 1    , & ! Order of the matrix
                  1                    , & ! Number of columns of the matrix R
                  DAL(ie_us + 1:ie_ds) , & ! The (n-1) sub-diagonal elements of A.
                  DA(ie_us:ie_ds)      , & ! The n diagonal elements of A
                  DAU(ie_us:ie_ds-1)   , & ! The (n-1) super-diagonal elements of A
                  R(ie_us:ie_ds)       , & ! right hand side
                  msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs + 1    , & ! Dimension of R
                  INFO )                                                    ! Info

      ! Update of the free surface height
      Hk(ie_us:ie_ds) = Hk(ie_us:ie_ds) + R(ie_us:ie_ds)
      !write(6,*) "test velocity",K, abs(dxHk)**0.5, Rh**(gam-1)
      !write(6,*) "test velocity", minval(v(ie_us+1:ie_ds-1)), maxval(v(ie_us+1:ie_ds-1))
      !write(6,*) "test height ", minval(Hk(ie_us:ie_ds)),maxval(Hk(ie_us:ie_ds))
      !pause
      ! update of the stopping criteria, the maximum of the height correction
      !write(6,*) "test output", maxval(abs(R(ie_us:ie_ds))),norm_err,it_sol
      norm_err  =  max( norm_err , sqrt(sum( R(ie_us:ie_ds)**2)) )

      ! update of the current iteration 
      it_sol=it_sol+1

   end do 
   !pause
   ! update of the height on ghost cell
   Hk( ie_us - 1 ) = Hk( ie_us ) + Hk( ie_us ) - Hk( ie_us + 1 )!+ msh%cs( ie_us )%bathy - msh%cs( ie_us + 1 )%bathy
   Hk( ie_ds + 1 ) = Hk( ie_ds ) + Hk( ie_ds ) - Hk( ie_ds - 1 ) !+ msh%cs( ie_ds )%bathy - msh%cs( ie_ds - 1 )%bathy
   do ie=ie_us ,ie_ds
      call perimeter(msh%cs(ie), Hk(ie) - msh%cs(ie)%bathy, pw)
      call strickler(msh%cs(ie), Hk(ie) - msh%cs(ie)%bathy, K)
      Rh = htoA(msh%cs(ie), Hk(ie) - msh%cs(ie)%bathy) / pw
      ! compute the velocity on the first and last cross-section
      if (ie == ie_us .or. ie == ie_ds) then 
         dxie =  2 * max(msh%cs(ie+1)%deltademi,msh%cs(ie)%deltademi)
      else 
         dxie =  ( msh%cs(ie+1)%deltademi+msh%cs(ie)%deltademi )        
      end if 
      dxHk = ( Hk(ie+1) - Hk(ie-1) ) / ( dxie )
      v(ie)=K * abs(dxHk)**0.5 * Rh**(gam-1)
      ! update of the unknown water depth
      dof%h(ie) = Hk(ie) - msh%cs(ie)%bathy
      ! update the unknown wet area and flow
      dof%A(ie)=htoA(msh%cs(ie), dof%h(ie))
      dof%Q(ie)= v(ie) * dof%A(ie)
      !write(6,*) "test output", Hk(ie+1),Hk(ie-1),dxie,abs(dxHk)**0.5,v(ie)
   end do  

     
END SUBROUTINE implicit_diffusive_wave
