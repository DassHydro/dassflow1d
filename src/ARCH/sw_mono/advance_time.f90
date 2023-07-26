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
!> \file advance_time.f90
!! \brief This file includes the computation of time step \f$dt\f$ from CFL and advancement of the time variable.

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Advancement of Time Variables / Eventual Adaptative Time Step Calculation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>  Advancement of time variables / Eventual adaptative time step calculation
!!
!! \details This subroutine compute the advancement of time variables and compute time step \f$dt\f$ from CFL value
!! if it is necessary.
!! \param[in]    dof Unknowns of the model.
!! \param[in]    msh Mesh of the problem.
SUBROUTINE advance_time( dof , msh )

   USE m_common
   USE m_mesh
   USE m_model
   USE m_time_screen                                                                                              !NOADJ

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE(Mesh), intent(in)  ::  msh

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp)  ::  c , vel , dist , s


!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   nt = nt + 1

   if ( nt == max_nt_for_direct ) end_time_loop = .true.

!======================================================================================================================!
!  Adaptative Time Step
!======================================================================================================================!
   if (adapt_dt ==1) then 

      dt=dt
!<NOADJ
      if ( fix_time_step_serie == 2 ) then
         read(80) tc , dt , end_time_loop
         return
      end if


      dt=hugem

      do i = 3,msh%ncs+2
         if ( dof%h(i) > 0.0_rp ) then

            dist = msh%cs(i)%delta
            c = sqrt( g * dof%h(i) )
            !c=sqrt(g*div_by_except_0(dof%s(i),surfaceWidth(i)))
            !c   =sqrt(g*div_by_except_0(dof%s(i),dof%w(i))
            vel = abs(div_by_except_0(dof%q(i),dof%s(i)))
            dt = min( dt , dist / ( vel + c ) )

         end if

      end do

      if (dt.eq.hugem) then
         dt=0.0001_rp
         tc=ts+100
      end if

      dt=cfl*dt

      if ( tc + dt < ts ) then
         tc = tc + dt

      else

         dt = ts - tc ; tc  = ts

         end_time_loop = .true.

      end if

      call write_scalar_in_time( dt , 'time_step' )

      if ( fix_time_step_serie == 1  ) then 
         write(80) tc , dt , end_time_loop
      endif
!>NOADJ
!======================================================================================================================!
!  Constant Time Step
!======================================================================================================================!
   else 
      tc = real(nt,rp) * dt

      if ( tc + dt + zerom > ts ) end_time_loop = .true.

   end if

!======================================================================================================================!
!  Output Time Informations
!======================================================================================================================!

   call Print_Screen ( 'dt') !NOADJ           
END SUBROUTINE advance_time
