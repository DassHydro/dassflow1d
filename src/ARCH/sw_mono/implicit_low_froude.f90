!======================================================================================================================!
!
!                    DassFlow Version 2.0
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
!> \file implicit_low_froude.f90
!! \brief This file includes the computation with implicit low froude methode without MUSCL and slope limiter
!! \details For more details, see developper documentation. This file contains only implicit_low_froude subroutine


!> Subroutine of the computation implicit low froude method without MUSCL and without slope limiter. 
!!
!! \details More details given in developper documentation
!! \param[in]  dof Unknowns of the model.
!! \param[in]    msh Mesh of the model.

#ifdef DEBLOCK
SUBROUTINE implicit_low_froude( dof , msh )

   USE m_common
   USE m_mesh
   USE m_time_screen
   USE m_model
   USE m_linear_solver
   USE m_user_data
   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE(Mesh), intent(in   )  ::  msh
   TYPE( unk ), intent(inout)  ::  dof

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp) :: norm_err    !Convergence criteria 
   integer(ip) :: it_sol   !Newton s method iteration
   integer(ip) :: n       ! Dimension of problem
   integer(ip) :: info,ncp2
   real(rp), dimension(msh%ncs-1) :: DL  ! Sub-diagonal elements of A
   real(rp), dimension(msh%ncs-1) :: DU  ! Sup-diagonal elements of A
   real(rp), dimension(msh%ncs) :: D   ! Diagonal elements of A

   real(rp), dimension(msh%ncs) :: B   ! Definition of B

   real(rp) :: wjm1,wj,wjp1
   real(rp) :: vjm1,vj,vjp1
   real(rp) :: hjm1,hj,hjp1
   real(rp) :: sjm1,sj,sjp1

   real(rp) :: GetWFromH

   real(rp), dimension( msh%ncs+4)  ::  h,v,s                    ! Temporal depths

   
   real(rp), dimension(msh%ncs-1) :: DL1  ! Sub-diagonal elements of A
   real(rp), dimension(2) :: DU1  ! Sup-diagonal elements of A
   real(rp), dimension(3) :: D1   ! Diagonal elements of A

   real(rp), dimension(3) :: B1   ! Definition of B


   real(rp) :: gamma_imp
   
   real(rp) :: Fjm1d2 , Fjp1d2,Rj

   real(rp) :: SourceTermPente,SourceTermFriction

   real(rp) ::A ,q,dz,HtoS,s_new,v_new,F2

   real(rp) :: manning, hydraulicradiusie

   !======================================================================================================================!
   !  Begin Subroutine
   !======================================================================================================================!

   call calc_boundary_state( msh,dof)
   !call SurfaceToHeight(msh,dof)

   !===================================================================================================================!
   !  Newton s method for nonlinear systems
   !===================================================================================================================!
   norm_err = hugem

   it_sol = 0
   
   h(:)  =  dof%h(:)
   v(:)  =  dof%q(:)/dof%s(:)
   s(:)  =  dof%s(:)

   h(2)=h(3)
   v(2)=dof%q(2)/dof%s(3)


   v(msh%ncs+3)=v(msh%ncs+2)
   
   
   gamma_imp=(2*g)/((msh%cs(10)%delta)*(msh%cs(10)%delta))

   do while ( it_sol <= 100 .and. norm_err > 10._rp * zerom )
      
      do ie= 4,msh%ncs+1

         !gamma_imp=(2*g)/((msh%cs(ie)%delta)*(msh%cs(ie)%delta))
         gamma_imp=(2*g)/((msh%cs(ie)%delta))

   
         vjm1  =v(ie-1) !v_(j-1)
         vj    =v(ie  ) !v_(j  )
         vjp1  =v(ie+1) !v_(j+1)

         hjm1  =dof%h(ie-1) !v_(j-1)
         hj    =dof%h(ie  ) !v_(j  )
         hjp1  =dof%h(ie+1) !v_(j+1)

         wjm1  =GetWFromH(msh,hjm1,ie-1)!msh%cs(ie-1)%surfacewidth(1)  !w_(j-1)
         wj    =GetWFromH(msh,hj  ,ie  )!msh%cs(ie  )%surfacewidth(1)  !w_(j  )
         wjp1  =GetWFromH(msh,hjp1,ie+1)!msh%cs(ie+1)%surfacewidth(1)  !w_(j+1)

         sjm1  =HtoS(msh,ie-1, hjm1) !hjm1*wjm1 !HtoS( msh,im, height )
         sj    =HtoS(msh,ie  , hj  ) !hj*wj
         sjp1  =HtoS(msh,ie+1, hjp1) !hjp1*wjp1

         ! Definition of DL
         DL(ie-3)=demi*(dt/msh%cs(ie)%delta)*(-vjm1*wjm1-demi*gamma_imp*dt*(sjm1+sj-(hj-hjm1)*wjm1))

         ! Definition of DU
         DU(ie-2)=demi*(dt/msh%cs(ie)%delta)*(vjp1*wjp1-demi*gamma_imp*dt*(sjp1+sj+(hjp1-hj)*wjp1))

         ! Defintion of D
         D(ie-2)=wj-d1p4*((gamma_imp*dt*dt)/msh%cs(ie)%delta)*(-(sj+sjp1)+(hjp1-hj)*wj-(sjm1+sj)-(hj-hjm1)*wj)

         ! Defintion of B
         B(ie-2)=sj-h(ie)*wj+demi*(dt/msh%cs(ie)%delta)*&
         ((vjp1*sjp1-vjm1*sjm1)-demi*gamma_imp*dt*((hjp1-hj)*(sj+sjp1)-(hj-hjm1)*(sjm1+sj)))


      end do


      !ie=3

      !gamma_imp=(2*g)/((msh%cs(3)%delta)*(msh%cs(3)%delta))      
      gamma_imp=(2*g)/((msh%cs(3)%delta))

      vjm1  =v(2) !v_(j-1)
      vj    =v(3) !v_(j  )
      vjp1  =v(4) !v_(j+1)

      hjm1  =dof%h(2) !v_(j-1)
      hj    =dof%h(3) !v_(j  )
      hjp1  =dof%h(4) !v_(j+1)

      wjm1  =GetWFromH(msh,hjm1,2)!msh%cs(ie-1)%surfacewidth(1)  !w_(j-1)
      wj    =GetWFromH(msh,hj  ,3)!msh%cs(ie  )%surfacewidth(1)  !w_(j  )
      wjp1  =GetWFromH(msh,hjp1,4)!msh%cs(ie+1)%surfacewidth(1)  !w_(j+1)

      sjm1  =HtoS(msh,2,hjm1) !hjm1*wjm1 !HtoS( msh,im, height )
      sj    =HtoS(msh,3,hj  ) !hj*wj
      sjp1  =HtoS(msh,4,hjp1) !hjp1*wjp1


      DU(1)=demi*(dt/msh%cs(ie)%delta)*(vjp1*wjp1-demi*gamma_imp*dt*(sjp1+sj+(hjp1-hj)*wjp1))
      D(1)=wj-d1p4*((gamma_imp*dt*dt)/msh%cs(ie)%delta)*(-(sj+sjp1)+(hjp1-hj)*wj-(sjm1+sj)-(hj-hjm1)*wj)
      B(1)=sj-h(3)*wj+demi*(dt/msh%cs(ie)%delta)*&
      ((vjp1*sjp1-vjm1*sjm1)-demi*gamma_imp*dt*((hjp1-hj)*(sj+sjp1)-(hj-hjm1)*(sjm1+sj)))


      
      !ie =msh%ncs+2
      ncp2=msh%ncs+2

      !gamma_imp=(2*g)/((msh%cs(ncp2)%delta)*(msh%cs(ncp2)%delta))
      gamma_imp=(2*g)/((msh%cs(ncp2)%delta))

      vjm1  =v(msh%ncs+1) !v_(j-1)
      vj    =v(msh%ncs+2) !v_(j  )
      vjp1  =v(msh%ncs+3) !v_(j+1)

      hjm1  =dof%h(msh%ncs+1) !v_(j-1)
      hj    =dof%h(msh%ncs+2) !v_(j  )
      hjp1  =dof%h(msh%ncs+3) !v_(j+1)

      wjm1  =GetWFromH(msh,hjm1,msh%ncs+1)!msh%cs(ie-1)%surfacewidth(1)  !w_(j-1)
      wj    =GetWFromH(msh,hj  ,msh%ncs+2)!msh%cs(ie  )%surfacewidth(1)  !w_(j  )
      wjp1  =GetWFromH(msh,hjp1,msh%ncs+3)!msh%cs(ie+1)%surfacewidth(1)  !w_(j+1)

      sjm1  =HtoS(msh,msh%ncs+1,hjm1) !hjm1*wjm1 !HtoS( msh,im, height )
      sj    =HtoS(msh,msh%ncs+2,hj  ) !hj*wj
      sjp1  =HtoS(msh,msh%ncs+3,hjp1) !hjp1*wjp1


      DU(msh%ncs-1)=demi*(dt/msh%cs(ie)%delta)*(-vjm1*wjm1-demi*gamma_imp*dt*(sjm1+sj-(hj-hjm1)*wjm1))
      D(msh%ncs)=wj-d1p4*((gamma_imp*dt*dt)/msh%cs(ie)%delta)*(-(sj+sjp1)+(hjp1-hj)*wj-(sjm1+sj)-(hj-hjm1)*wj)
      B(msh%ncs)=sj-h(msh%ncs+2)*wj+demi*(dt/msh%cs(ie)%delta)*&
      ((vjp1*sjp1-vjm1*sjm1)-demi*gamma_imp*dt*((hjp1-hj)*(sj+sjp1)-(hj-hjm1)*(sjm1+sj)))


      norm_err  =  maxval( abs( B(:) ) )

      !Call of linear solver (Solve AX=B with A tridiagonal matrix)
      call DGTSV( msh%ncs    , & ! Order of the matrix
                  1           , & ! Number of columns of the matrix B
                  DL          , & ! The (n-1) sub-diagonal elements of A.
                  D           , & ! The n diagonal elements of A
                  DU          , & ! The (n-1) super-diagonal elements of A
                  B           , & ! B
                  msh%ncs    , & ! Dimension of B
                  INFO )          ! Info
     
      dof%h(3:msh%ncs+2)  =  dof%h(3:msh%ncs+2)  - B(:)

      !Boundary conditions
      dof%h(2)=dof%h(3)
      dof%s(2)=dof%s(3)
      
      norm_err  =  min( norm_err , maxval( abs( B(:) ) ) )

      it_sol=it_sol+1
      
   end do

   if (it_sol.ge.100) then
      write(6,'(I3,ES15.7)') it_sol , norm_err
   end if  
   


   !HtoS
   do ie=1,msh%ncs+4
      dof%s(ie)=HtoS( msh,ie, dof%h(ie) )
   end do


   dof%s(2)=dof%s(3)
   s(2)=s(3)

   !===================================================================================================================!
   !  Calculating explicit Q
   !===================================================================================================================!
   do ie=3,msh%ncs+2

      vjm1=v( ie-1 )
      vj  =v( ie   )
      vjp1=v( ie+1 )
   
      Fjm1d2=demi*(v(ie-1)*dof%s(ie-1)+v(ie)*dof%s(ie))-d1p4*gamma_imp*dt*(dof%h(ie  )-dof%h(ie-1))*(dof%s(ie-1)+dof%s(ie))
      Fjp1d2=demi*(v(ie+1)*dof%s(ie+1)+v(ie)*dof%s(ie))-d1p4*gamma_imp*dt*(dof%h(ie+1)-dof%h(ie  ))*(dof%s(ie+1)+dof%s(ie))

      dz = (bathy_cell(ie+1)-bathy_cell(ie-1))&
            /(msh%cs(ie)%deltademi+msh%cs(ie-1)%deltademi) !Slope 

      SourceTermPente =-g*dz*dof%s(ie)!*wj                        !Slope          term

      q=dof%q(ie)

      !===================================================
      !  S,Q computation with source term (friction)
      !===================================================
      if (friction.eq.1) then 
         call computationRhCrossSection(msh,dof,hydraulicRadiusie,ie)        !Hydraulics computation (if friction==1)
            call calc_K_at_cs(dof, Manning, ie)                  ! Update of Manning
         A=(q*abs(q))/(((Manning**2))*(dof%s(ie))*(hydraulicRadiusie**d4p3))  ! Sf= |q|*q /(K^2*S*Rh^(4/3)                
         SourceTermFriction=-g*A
      else 
         SourceTermFriction=0._rp
      end if


      !Rj      = SourceTermPente+SourceTermFriction   !Pressure+slope term+Friction Term
      Rj      = SourceTermFriction
      !Rj      = 0.0_rp

      v_new=vj*s(ie)-(dt/msh%cs(ie)%delta)*&
      ((max(Fjp1d2,0._rp)*vj+min(Fjp1d2,0._rp)*vjp1)-(max(Fjm1d2,0._rp)*vjm1+min(Fjm1d2,0._rp)*vj)) &
      !((max(Fjp1d2,0._rp)*vjp1+min(Fjp1d2,0._rp)*vj)-(max(Fjm1d2,0._rp)*vj+min(Fjm1d2,0._rp)*vjm1)) &
      -((g*dt)/(2._rp*msh%cs(ie)%delta))*dof%s(ie)*&
      ((dof%h(ie+1)+bathy_cell(ie+1))-(dof%h(ie-1)+bathy_cell(ie-1)))&
      +dt*Rj


      if ((v_new.le.zerom)) then
         !write(*,*) dof%q(ie)
         v_new=0.0
      endif

      dof%q(ie)=v_new !*dof%s(ie)

   end do


END SUBROUTINE implicit_low_froude

#endif