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
!> \file explicit_low_froude_b2.f90
!! \brief This file includes the computation with explicit low froude method with MUSCL and with slope limiter.
!! An unique subroutine :
!! explicit_low_froude

!> Subroutine of the computation explicit low froude method with MUSCL and with slope limiter. 
!!
!! \details Algorithm used:
!!
!!     Update boundaries conditions
!!     Update mvector
!!
!!     do i = 2, n:
!!       Computation of MUSCL variable ( h and u) with slope limiter
!!       Computation of F (low froude)
!!       s(i) = max(0,s-(dt/dxi) (F1(i+1)-F1(i)))
!!       Update mvector(i)
!!       Update H(i) (from S(i))
!!     end do
!!     
!!     do i = 2, n:
!!       Computation of Sg(i) and Sf(i)
!!       Computation of q(i)
!!     end do
!!
!! \param[in]  dof Unknowns of the model.
!! \param[in]    msh Mesh of the model.SUBROUTINE explicit_low_froude_b2( dof , msh )

#ifdef DEBLOCK
SUBROUTINE explicit_low_froude_b2( dof , msh )
   USE m_common
   USE m_mesh
   USE m_time_screen          !NOADJ
   USE m_model
   USE m_linear_solver        !NOADJ
   USE m_user_data
   USE m_numeric
   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE(Mesh), intent(in   )  ::  msh
   TYPE( unk ), intent(inout)  ::  dof

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp) :: wjm1,wj,wjp1
   real(rp) :: vjm1,vj,vjp1
   real(rp) :: hjm1,hj,hjp1
   real(rp) :: sjm1,sj,sjp1

   real(rp), dimension( msh%ncs+4)  ::  h,v,s                    ! Temporal depths

   real(rp) :: gamma_expl,alpha_expl
   
   real(rp) :: Fjm1d2 , Fjp1d2,Phijp1d2, Phijm1d2,Rj

   real(rp) :: SourceTermPente,SourceTermFriction

   real(rp) ::A ,q,dz,HtoS,s_new,v_new,F2

   real(rp) :: manning, hydraulicradiusie



   real(rp)  :: hm2,hm1,hi,hp1,hp2   !Height    at the point i-2,i-1,i,i+1,i+2

   real(rp)  :: qm2,qm1,qi,qp1,qp2   !Flow      at the point i-2,i-1,i,i+1,i+2

   real(rp)  :: sm2,sm1,si,sp1,sp2   !Surface   at the point i-2,i-1,i,i+1,i+2

   real(rp)  :: um2,um1,ui,up1,up2   !Velocity  at the point i-2,i-1,i,i+1,i+2

   real(rp)  :: zm2,zm1,zi,zp1,zp2   !Elevation at the point i-2,i-1,i,i+1,i+2


   integer(ip), dimension(msh%ncs+4) :: mVector              ! Pressure for fluxes computing

   real(rp) :: xm1demiMxm1,xm1demiMxm3demi,xm1demiMx,xp1demiMxm1demi                   ! Deltax value  
   
   real(rp) ::deltaxm1demi, deltax,deltaxm1                                            !  Deltax and Deltax1demi

   real(rp) :: deltahath,deltahathm1,sigmah,sigmahm1,deltah,deltahm1,hm1demim,hm1demip ! Delta and sigma h for h reconstruction

   real(rp) :: deltahatu,deltahatum1,sigmau,sigmaum1,deltau,deltaum1,um1demim,um1demip ! Delta and sigma u for u reconstruction

   real(rp) :: deltahats,deltahatsm1,sigmas,sigmasm1,deltas,deltasm1,sm1demim,sm1demip ! Delta and sigma u for u reconstruction

   real(rp)  :: sl,sr ! Surface  to the cell left and right

   real(rp)  :: ul,ur ! Velocity to the cell left and right

   real(rp)  :: ql,qr ! Flow     to the cell left and right

   real(rp)  :: hl,hr ! Height   to the cell left and right

   !real(rp), dimension( 2,msh%ncs+4 )  ::  tflux    ! Vector total flux (aggregation of flux)
   real(rp), dimension(msh%ncs+4 )  ::   tphi

   real(rp) :: L

   real(rp) :: HtoSInterface                          ! Transformation of H to S to the interface between two cell

   real(rp) :: alphab                                 ! Alpha (For MUSCL Reconstruction)

   real(rp) :: zm1demi                                ! Value of z to the interface between cell i-1 and i

   real(rp) :: deltaxm,deltaxp                        ! DeltaX at cell i-1, DeltaX at cell i+1


   real(rp) :: GetWFromH

   !======================================================================================================================!
   !  Begin Subroutine
   !======================================================================================================================!
   !write(*,*) 'expli'
   call calc_boundary_state( msh,dof)

   call UpdateMVector(msh,dof,mvector)         ! Update mvector   

   call SurfaceToHeightCrossSection(msh,dof,1,mvector)          ! Update of dof%h(1)
   call SurfaceToHeightCrossSection(msh,dof,2,mvector)          ! Update of dof%h(2)
   call SurfaceToHeightCrossSection(msh,dof,msh%ncs+3,mvector) ! Update of dof%h(msh%ncs+3)
   call SurfaceToHeightCrossSection(msh,dof,msh%ncs+4,mvector) ! Update of dof%h(msh%ncs+4)


      
   h(:)  =  dof%h(:)
   do ie=1,msh%ncs+4
      v(ie)  =  dof%q(ie)/dof%s(ie)
   end do 
   s(:)  =  dof%s(:)

!    h(1)=h(3)
!    h(2)=h(3)
!    v(1)=dof%q(2)/dof%s(3)
!    v(2)=dof%q(2)/dof%s(3)


   v(msh%ncs+3)=v(msh%ncs+2)
   v(msh%ncs+4)=v(msh%ncs+2)


   ! Original 
   gamma_expl=0.5_rp
   alpha_expl=1._rp 
   L=1._rp 

   !For Garonne :   
   !gamma_expl=1._rp
   !alpha_expl=2._rp 
   !L=5._rp 


   gamma_expl=0.1_rp
   alpha_expl=5._rp 
   L=1000._rp 


   !write(*,*) 'dde'
   do ie=3,msh%ncs+3

     !===================================
      !  variable definition  
      !===================================
      !surfaceWidth(ie)=getWFromH(msh,dof%h(ie),ie)

      zm2 =bathy_cell(ie-2)         
      zm1 =bathy_cell(ie-1)
      zi  =bathy_cell(ie  )
      zp1 =bathy_cell(ie+1)
!       if (tc == tc0 + dt) print '(I4,4E12.5)', ie, zm2, zm1, zi, zp1

      hm2 =dof%h(ie-2)
      hm1 =dof%h(ie-1)
      hi   =dof%h(ie  )
      hp1 =dof%h(ie+1)

      hm2=hm2+zm2
      hm1=hm1+zm1
      hi  =hi  +zi
      hp1=hp1+zp1

      qm2 =dof%q(ie-2)
      qm1 =dof%q(ie-1)
      qi   =dof%q(ie  )
      qp1 =dof%q(ie+1)

      sm2 =dof%s(ie-2)      
      sm1 =dof%s(ie-1)
      si  =dof%s(ie  )
      sp1 =dof%s(ie+1)

      deltaxm1     = msh%cs(ie-1)%delta
      deltax       = msh%cs(ie)%delta
      deltaxm1demi = msh%cs(ie)%deltademi

      um2 =div_by_except_0(qm2,sm2)
      um1 =div_by_except_0(qm1,sm1)
      ui  =div_by_except_0(qi ,si )
      up1 =div_by_except_0(qp1,sp1)



      !===================================
      !  Interface value + MUSCL (see doc)
      !===================================
      alphab=one
      xm1demiMxm1=demi*deltaxm1demi
      xm1demiMxm3demi=deltaxm1

      xm1demiMx=-xm1demiMxm1
      xp1demiMxm1demi=deltax

      !  deltahat computing
      deltahath=demi*(hp1-hm1)
      deltahathm1=demi*(hi-hm2)

      deltahatu=demi*(up1-um1)
      deltahatum1=demi*(ui-um2)

      deltahats=demi*(sp1-sm1)
      deltahatsm1=demi*(si-sm2)


      !  sigma computing
      sigmah  =div_by_except_0(deltahath  ,abs(deltahath  ))
      sigmau  =div_by_except_0(deltahatu  ,abs(deltahatu  ))
      sigmas  =div_by_except_0(deltahats  ,abs(deltahats  ))

      sigmahm1=div_by_except_0(deltahathm1,abs(deltahathm1))
      sigmaum1=div_by_except_0(deltahatum1,abs(deltahatum1))
      sigmasm1=div_by_except_0(deltahatsm1,abs(deltahatsm1))

      
      !  deltah computing
      deltah  =sigmah  *max(zero, min(sigmah  *alphab*(hi -hm1),abs(deltahath  ),sigmah  *alphab*(hp1-hi )))  
      deltahm1=sigmahm1*max(zero, min(sigmahm1*alphab*(hm1-hm2),abs(deltahathm1),sigmahm1*alphab*(hi -hm1)))  
      deltau  =sigmau  *max(zero, min(sigmau  *alphab*(ui -um1),abs(deltahatu  ),sigmau  *alphab*(up1-ui )))  
      deltaum1=sigmaum1*max(zero, min(sigmaum1*alphab*(um1-um2),abs(deltahatum1),sigmaum1*alphab*(ui -um1)))  
    
      deltas  =sigmas  *max(zero, min(sigmas  *alphab*(si  -sm1),abs(deltahats ),sigmas  *alphab*(sp1-si  )))  
      deltasm1=sigmasm1*max(zero, min(sigmasm1*alphab*(sm1-sm2),abs(deltahatsm1),sigmasm1*alphab*(si  -sm1)))  

      !  hm1demim & hm1demip
      deltaxm=(xm1demiMxm1)/(xm1demiMxm3demi)
      deltaxp=(xm1demiMx  )/(xp1demiMxm1demi)

      hm1demim=hm1 + deltahm1*deltaxm
      hm1demip=hi  + deltah *deltaxp

      um1demim=um1 + deltaum1*deltaxm
      um1demip=ui  + deltau *deltaxp

      sm1demim=sm1 + deltasm1*deltaxm
      sm1demip=si  + deltas  *deltaxp

      !===================================
      !  transformation (h,u) to (S,Q) 
      !===================================
      Sl=HtoSInterface(msh,ie-1,ie,hm1demim)
      Sr=HtoSInterface(msh,ie-1,ie,hm1demip)
      !sl=sm1demim
      !sr=sm1demip
            

      !zm1demi=demi*(zm1+zi)
      !hm1demim=max(zero,hm1demim-zm1demi) 
      !hm1demip=max(zero,hm1demip-zm1demi)

         
      Ql=sl*um1demim
      Qr=sr*um1demip


      !L=1.0_rp !5.0_rp 
      !write(*,*) L

      tflux1(ie)=demi*(um1demim*sl+um1demip*sr)-d1p4*gamma_expl*&
      (dt/msh%cs(ie)%delta)*g*(hm1demip-hm1demim)*(Sl+Sr)

      if (ie.eq.3) then
         tflux1(ie)=dof%q(1) !linear_interp( bc%hyd%t ,bc%hyd%q ,tc) !65.0_rp
      end if
      !write(*,*) 'ddd'
      zm1demi=demi*(zm1+zi)
      hm1demim=max(zero,hm1demim-zm1demi) 
      hm1demip=max(zero,hm1demip-zm1demi)

!       tphi(ie)=demi*(hm1+hi)-(alpha_expl/L)*(dt/msh%cs(ie)%delta)*demi*(um1demip*hm1demip-um1demim*hm1demim)
!       !tphi(ie) = demi*(h(ie-1)+bathy_cell(ie-1)+h(ie)+bathy_cell(ie))-(alpha_expl/L)*&
!       !(dt/msh%cs(ie)%delta)*demi*(v(ie  )*(s(ie  ))-v(ie-1)*(s(ie-1)))
      
      tphi(ie)=demi*(hm1+hi)-(alpha_expl/(2*L))*(dt/msh%cs(ie)%delta)*demi*(um1demip*sr-um1demim*sl)
      
      
   end do
!    if (tc == tc0 + dt) read(*,*)


   !tflux1(3) = 65.0_rp
   !===================================================================================================================!
   !  Calculating explicit Q
   !===================================================================================================================! 
   do ie=3,msh%ncs+2
      vjm1=v( ie-1 )
      vj  =v( ie   )
      vjp1=v( ie+1 )
   
      Fjm1d2=tflux1(ie)
      Fjp1d2=tflux1(ie+1)

      dof%s(ie)=dof%s(ie)-(dt/msh%cs(ie)%delta)*(Fjp1d2-Fjm1d2)

      Phijm1d2=tphi(ie)
      Phijp1d2=tphi(ie+1)

      q=dof%q(ie)

      call UpdateMVectorElement(msh,dof,ie,mvector)                 ! Update of mvector 
      !!write(*,*) 'coucou'
      call SurfaceToHeightCrossSection(msh,dof,ie,mvector)          ! Update of dof%h(ie)

      !===================================================
      !  S,Q computation with source term (friction)
      !===================================================
      if (friction.eq.1) then 
         call computationRhCrossSection(msh,dof,hydraulicRadiusie,ie,mvector) !Hydraulics computation (if friction==1)
            call calc_K_at_cs(dof, Manning, ie)                  ! Update of Manning
         A=(q*abs(q))/(((Manning**2))*(s(ie))*(hydraulicRadiusie**d4p3))  ! Sf= |q|*q /(K^2*S*Rh^(4/3)                
         SourceTermFriction=-g*A
      else 
         SourceTermFriction=0._rp
      end if

      Rj      = SourceTermFriction
      
      dof%q(ie)=v(ie)*s(ie)-(dt/msh%cs(ie)%delta)*&
      ((max(Fjp1d2,0._rp)*vj+min(Fjp1d2,0._rp)*vjp1)-(max(Fjm1d2,0._rp)*vjm1+min(Fjm1d2,0._rp)*vj)) &
      -((g*dt)/(msh%cs(ie)%delta))*s(ie)*(Phijp1d2-Phijm1d2)&
      +dt*Rj
      if ((abs(dof%q(ie)).le.(0.00000000000001))) then
         dof%q(ie)=0.0
      endif
      !dof%q(ie)=65.0_rp



   end do


END SUBROUTINE explicit_low_froude_b2
#endif
