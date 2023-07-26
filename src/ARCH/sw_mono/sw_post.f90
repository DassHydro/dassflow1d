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
!> \file sw_post.f90
!! \brief This file include post-processing subroutine.
!! Only one subroutine is in this file : sw_post_treatment subroutine.

!> This subroutine done the post-processing of results.
!! \param[in] msh  Mesh of the model.
!! \param[in] dof Unknowns of the model.
SUBROUTINE sw_post_treatment( dof , msh )

   USE m_common
   USE m_mesh
   USE m_time_screen                                                                                              !NOADJ
   USE m_model
   USE m_user_data

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type(Mesh), intent(in)  ::  msh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!
                                                                                                                 !<NOADJ
   real(rp)  ::  mass

   real(rp)  ::  norm_L2_e(2)
   real(rp)  ::  norm_L2_s(2)

   real(rp)  ::  s_e , q_e , diff , surf_total

   real(rp)  ::  discharg , mass_cut_tot                                                                         !>NOADJ

   integer(ip)  ::  num_bc

                                                                                                                 !<NOADJ

   call write_stations( dof, msh )
!======================================================================================================================!
!  Norms infinity, L1 and L2 in case of validation Makefile flag and a prescribed exact solution
!======================================================================================================================!

#ifdef USE_VALID

      !================================================================================================================!
      !  Output Norms in time or integrated in time or at the end of the simulation (t = ts)
      !================================================================================================================!

      norm_L2_s(:)=norm_L2 (:) 
      norm_L2 (:)  =  0._rp ; norm_L2_e (:)  =  0._rp;

      
      surf_total  =  0._rp

      do i = 1,msh%ncs

         !==========================================================================================================!

         s_e  =  max( 0._rp , s_exact( msh%crossSection(i)%coord%x , msh%crossSection(i)%coord%y , tc ) )
         q_e  =  q_exact( msh%crossSection(i)%coord%x , msh%crossSection(i)%coord%y , tc )

         !==========================================================================================================!

         diff  =  abs( dof%s(i) - s_e )

         norm_L2   (1)  =  norm_L2  (1)  +  msh%crossSection(i)%delta * diff**2
         norm_L2_e (1)  =  norm_L2_e(1)  +  msh%crossSection(i)%delta * s_e**2

         !==========================================================================================================!

         diff  =  abs( dof%q(i) - q_e )

         norm_L2   (2)  =  norm_L2  (2)  +  msh%crossSection(i)%delta * diff**2
         norm_L2_e (2)  =  norm_L2_e(2)  +  msh%crossSection(i)%delta * q_e**2

         !==========================================================================================================!

         !surf_total  =  surf_total  +  msh%crossSection(i)%delta
      end do


      norm_L2 (1)  =  sqrt( div_by_except_0( norm_L2 (1) , norm_L2_e (1) ) )
      norm_L2 (2)  =  sqrt( div_by_except_0( norm_L2 (2) , norm_L2_e (2) ) )


      norm_L2 (1)  =  norm_L2 (1)  +  norm_L2_s (1)
      norm_L2 (2)  =  norm_L2 (2)  +  norm_L2_s (2)
      if ( end_time_loop ) then
               norm_L2 (:)  =  norm_L2 (:) / real( nt , 8 )
      end if

      call write_scalar_in_time( norm_L2 (1) , 'norm_rel_s_L2'  )
      call write_scalar_in_time( norm_L2 (2) , 'norm_rel_q_L2'  )


      if (end_time_loop ) then
         call write_pscalar( (/ norm_L2(1) ,norm_L2(2) /),&
                             (/ 'norm_rel_L2_s' ,'norm_rel_L2_q' /),&
                            'norms_rel_int_time' , 2 )
      end if

      if ( end_time_loop ) then

            call Print_Screen( 'norms_s' )
            call Print_Screen( 'norms_q' )
!            call Print_Screen( 'norms_u' )

      end if


#endif                                                                                                   
!>NOADJ

END SUBROUTINE sw_post_treatment
