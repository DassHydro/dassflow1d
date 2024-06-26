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
!> \file m_tap_vars.f90
!! \brief This file includes m_tap_vars module.
!! \details The file includes only m_tap_vars module (see doc m_tap_vars module).


!> Module m_tap_vars.
!!
!! \details This module allocates necessary variables to run linear tangent model and adjoint model of DassFlow
MODULE m_tap_vars

   USE m_common
   USE m_obs

   #ifdef USE_SW_MONO
      USE m_model
   #endif

   real(rp)  ::  dt_diff         !< Declaration of the derivative of \f$dt\f$.
   real(rp)  ::  dt_back         !< Declaration of the adjoint of \f$dt\f$.
   real(rp)  ::  tc_diff         !< Declaration of the derivative of \f$T\f$.
   real(rp)  ::  tc_back         !< Declaration of the adjoint of \f$T\f$.

   type( bcs )  ::  bc_diff      !< Derivative of boundaries conditions.
   type( bcs )  ::  bc_back      !< Adjoint of boundaries conditions.

   type( innovation_obs ), dimension(:), allocatable  ::  innovation_diff !< Derivative of innonvation vector (cost function).
   type( innovation_obs ), dimension(:), allocatable  ::  innovation_back !< Adjoint of innonvation vector (cost function).




   real(rp), dimension(:), allocatable  ::  tflux1_diff !< Derivative of flux1 vector (used to accelerate HLL).
   real(rp), dimension(:), allocatable  ::  tflux1_back !< Adjoint of flux2 vector (used to accelerate HLL).

   real(rp), dimension(:), allocatable  ::  tflux2_diff !< Derivative of flux2 vector (used to accelerate HLL).
   real(rp), dimension(:), allocatable  ::  tflux2_back !< Adjoint of flux2 vector (used to accelerate HLL).

   #ifdef USE_SW_MONO

   real(rp) :: alpha_ratcurve_diff !< Derivative of \f$\alpha_r\f$ (ratcurve model \f$h=\alpha_r Q ^{\beta_r}\f$).
   real(rp) :: alpha_ratcurve_back !< Adjoint of \f$\alpha_r\f$ (ratcurve model \f$h=\alpha_r Q ^{\beta_r}\f$).
   real(rp) :: beta_ratcurve_diff  !< Derivative of \f$\beta_r\f$ (ratcurve model \f$h=\alpha_r Q ^{\beta_r}\f$).
   real(rp) :: beta_ratcurve_back  !< Adjoint of \f$\beta_r\f$ (ratcurve model \f$h=\alpha_r Q ^{\beta_r}\f$).


   real(rp), dimension(:), allocatable  ::   bathy_cell_diff !< Derivative of \f$z_b\f$.
   real(rp), dimension(:), allocatable  ::   bathy_cell_back !< Adjoint of \f$z_b\f$.
   real(rp), dimension(:), allocatable  ::   bathy_points_diff
   real(rp), dimension(:), allocatable  ::   bathy_points_back
   type(strickler_params) :: K_params_diff                   !< Derivative of the Strickler parameters
   type(strickler_params) :: K_params_back                   !< Adjoint of the Strickler parameters

   real(rp), dimension(:), allocatable  ::   bathy_cell_chg_diff 
   real(rp), dimension(:), allocatable  ::   bathy_cell_chg_back 
   real(rp), dimension(:), allocatable  ::   bathy_points_chg_diff
   real(rp), dimension(:), allocatable  ::   bathy_points_chg_back
   
   real(rp), dimension(:), allocatable  ::   alpha_K_chg_diff 
   real(rp), dimension(:), allocatable  ::   alpha_K_chg_back 
   real(rp), dimension(:), allocatable  ::   beta_K_chg_diff 
   real(rp), dimension(:), allocatable  ::   beta_K_chg_back

   real(rp), dimension(:), allocatable  ::   qin_chg_diff 
   real(rp), dimension(:), allocatable  ::   qin_chg_back 
   real(rp) :: alpha_regul_diff
   real(rp) :: alpha_regul_back
   
   real(rp), dimension(:,:), allocatable  ::   qlat_chg_diff 
   real(rp), dimension(:,:), allocatable  ::   qlat_chg_back

   #endif


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of needed variables to run the linear tangent model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   !> Allocate space for derivative variable (necessary to run the linear tangent model). 
   !!
   !! \details The subroutine allocate space for all derivative variable used in this model (BC, dof, ...).
   !! \param[inout] dof0_diff Derivative of initiales conditions.  
   !! \param[inout] dof_diff Derivative of unknowns of the model.
   !! \param[in] mesh Mesh of the model.
   SUBROUTINE alloc_diff_vars( dof0_diff , dof_diff , mesh )

      type( unk ), intent(inout)  ::  dof0_diff , dof_diff

      type(Mesh), intent(in)  ::  mesh

      ! dof
      call alloc_dof( dof0_diff , mesh )
      call alloc_dof( dof_diff  , mesh )

      ! innovation
      allocate( innovation_diff( size( innovation ) ) )

      do iobs = 1,size( innovation )
         allocate( innovation_diff( iobs )%diff( size( innovation( iobs )%diff ) ) )
         innovation_diff( iobs )%diff(:)  =  0._rp
      end do

      ! Hydro file_1 
      allocate( bc_diff%hyd%t( size( bc%hyd%t ) ) )
      allocate( bc_diff%hyd%q( size( bc%hyd%q ) ) )

      bc_diff%hyd%t(:)  =  0._rp
      bc_diff%hyd%q(:)  =  0._rp
     
      ! Hydro Fourier series
      bc_diff%hyd_FS%a0    =  0._rp

      allocate( bc_diff%hyd_FS%n( size( bc%hyd_FS%n ) ) )
      allocate( bc_diff%hyd_FS%a( size( bc%hyd_FS%a ) ) )
      allocate( bc_diff%hyd_FS%b( size( bc%hyd_FS%b ) ) )

      bc_diff%hyd_FS%n(:)  =  0_ip
      bc_diff%hyd_FS%a(:)  =  0._rp
      bc_diff%hyd_FS%b(:)  =  0._rp
      
      
      !Hyd lat
      allocate( bc_diff%hyd_lat%t( size( bc%hyd_lat%t )))
      allocate( bc_diff%hyd_lat%q( size( bc%hyd_lat%x ), size( bc%hyd_lat%t )))
      allocate( dof_diff%qlat( 2, size( bc%hyd_lat%x )))
      
      bc_diff%hyd_lat%t(:) = 0._rp
      bc_diff%hyd_lat%q(:,:) = 0._rp
      dof_diff%qlat(:,:) = 0._rp
      
      ! Downstream elevations
      allocate( bc_diff%dsz%t( size( bc%dsz%t ) ) )
      allocate( bc_diff%dsz%y( size( bc%dsz%y ) ) )

      bc_diff%dsz%t(:)  =  0._rp
      bc_diff%dsz%y(:)  =  0._rp

      ! ratcurve
      allocate( bc_diff%rat%h( size( bc%rat%h ) ) )
      allocate( bc_diff%rat%q( size( bc%rat%q ) ) )

      bc_diff%rat%h(:)  =  0._rp
      bc_diff%rat%q(:)  =  0._rp

      alpha_ratcurve_diff=0.0_rp
      beta_ratcurve_diff=0.0_rp

      !flux
      allocate( tflux1_diff(mesh%ncs+4) )
      allocate( tflux2_diff(mesh%ncs+4) )
      tflux1_diff(:)=0._rp
      tflux2_diff(:)=0._rp


      dt_diff = 0._rp

      #ifdef USE_SW_MONO
         ! bathy
         allocate( bathy_cell_diff( size( bathy_cell ) ) )
         bathy_cell_diff(:)  =  0._rp
         if (allocated(bathy_points)) then
            allocate( bathy_points_diff( size( bathy_points ) ) )
         end if 

         ! Strickler
         select case ( K_model )
         case ( 'K(h)' )
            allocate( K_params_diff%alpha(size(K_params%alpha)) )
            allocate( K_params_diff%beta(size(K_params%alpha)) )
            K_params_diff%alpha(:) = 0._rp
            K_params_diff%beta(:) = 0._rp
         case default
            ! do nothing
         end select
         
         ! derivatives of changed variables     
         allocate( bathy_cell_chg_diff(nb_bathy_control_pts) )
         bathy_cell_chg_diff(:)  =  0._rp
         if (allocated(bathy_points_chg_diff)) then
            allocate( bathy_points_diff( size( bathy_points ) ) )
            bathy_points_diff(:)  =  0._rp
         end if 
         
            
         select case ( K_model )
         case ( 'K(h)' )
            allocate( alpha_K_chg_diff(nb_diff_K_reachs) )
            allocate( beta_K_chg_diff(nb_diff_K_reachs) )
            alpha_K_chg_diff(:) = 0._rp
            beta_K_chg_diff(:) = 0._rp
         case default
            ! do nothing
            print*, 'K model not taken into account'
            stop
         end select
            
         allocate( qin_chg_diff(size( qin_chg )) )
         allocate( qlat_chg_diff(size( qlat_chg(:,1)), size(qlat_chg(1,:)) ))
   
#endif

   END SUBROUTINE alloc_diff_vars


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of needed variables to run the adjoint model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Allocate space for ajdoint variable (necessary to run the adjoint model). 
   !!
   !! \details The subroutine allocate space for all derivative variable used in this model (BC, dof, ...).
   !! \param[inout] dof0_back Adjoint of initiales conditions.  
   !! \param[inout] dof_back Adjoint of unknowns of the model.
   !! \param[in] mesh Mesh of the model.
   SUBROUTINE alloc_back_vars( dof0_back , dof_back , mesh )

      type( unk ), intent(inout)  ::  dof0_back , dof_back

      type(Mesh), intent(in)  ::  mesh

      ! dof
      call alloc_dof( dof0_back , mesh )
      call alloc_dof( dof_back  , mesh )

      ! innovation
      allocate( innovation_back( size( innovation ) ) )

      do iobs = 1,size( innovation )
         allocate( innovation_back( iobs )%diff( size( innovation( iobs )%diff ) ) )
         innovation_back( iobs )%diff(:)  =  0._rp
      end do

      ! Hyd file 
      allocate( bc_back%hyd%t( size( bc%hyd%t ) ) )
      allocate( bc_back%hyd%q( size( bc%hyd%q ) ) )

      bc_back%hyd%t(:)  =  0._rp
      bc_back%hyd%q(:)  =  0._rp

      ! Hyd Fourier serie
      bc_back%hyd_FS%a0    =  0._rp

      allocate( bc_back%hyd_FS%n( size( bc%hyd_FS%n ) ) )
      allocate( bc_back%hyd_FS%a( size( bc%hyd_FS%a ) ) )
      allocate( bc_back%hyd_FS%b( size( bc%hyd_FS%b ) ) )

      bc_back%hyd_FS%n(:)  =  0_ip
      bc_back%hyd_FS%a(:)  =  0._rp
      bc_back%hyd_FS%b(:)  =  0._rp

      !Hyd lat
      allocate( bc_back%hyd_lat%t( size( bc%hyd_lat%t )))
      allocate( bc_back%hyd_lat%q( size( bc%hyd_lat%x ), size( bc%hyd_lat%t )))
      allocate( dof_back%qlat( 2, size( bc%hyd_lat%x )))
      
      bc_back%hyd_lat%t(:) = 0._rp
      bc_back%hyd_lat%q(:,:) = 0._rp
      dof_back%qlat(:,:) = 0._rp
      
      ! Downstream elevations
      allocate( bc_back%dsz%t( size( bc%dsz%t ) ) )
      allocate( bc_back%dsz%y( size( bc%dsz%y ) ) )

      bc_back%dsz%t(:)  =  0._rp
      bc_back%dsz%y(:)  =  0._rp
      
      ! Rating curve
      allocate( bc_back%rat%h( size( bc%rat%h ) ) )
      allocate( bc_back%rat%q( size( bc%rat%q ) ) )

      bc_back%rat%h(:)  =  0._rp
      bc_back%rat%q(:)  =  0._rp

      alpha_ratcurve_back  = 0.0_rp
      beta_ratcurve_back   = 0.0_rp

      ! flux
      allocate( tflux1_back(mesh%ncs+4))
      allocate( tflux2_back(mesh%ncs+4))
      tflux1_back(:)=0._rp
      tflux2_back(:)=0._rp


      dt_back = 0._rp

#ifdef USE_SW_MONO

      ! bathy
      allocate( bathy_cell_back( size( bathy_cell ) ) )
      bathy_cell_back(:)  =  0._rp
      if (allocated(bathy_points)) then
        allocate( bathy_points_back( size( bathy_points ) ) )
      end if 

      ! Strickler
      select case ( K_model )
      case ( 'K(h)' )
         allocate( K_params_back%alpha(size(K_params%alpha)) )
         allocate( K_params_back%beta(size(K_params%alpha)) )
         K_params_back%alpha(:) = 0._rp
         K_params_back%beta(:) = 0._rp
      case default
         ! do nothing
      end select

      ! adjoint of changed variables            
      allocate( bathy_cell_chg_back(nb_bathy_control_pts) )
      bathy_cell_chg_back(:)  =  0._rp
      if (allocated(bathy_points)) then
        allocate( bathy_points_chg_back( size( bathy_points ) ) )
        bathy_points_chg_back(:)  =  0._rp
      end if 
      
      select case ( K_model )
      case ( 'K(h)' )
         allocate( alpha_K_chg_back(K_params%nb_diff_K_reachs) )
         allocate( beta_K_chg_back(K_params%nb_diff_K_reachs) )
         alpha_K_chg_back(:) = 0._rp
         beta_K_chg_back(:) = 0._rp
      case default
         ! do nothing
         print*, 'K model not taken into account'
         stop
      end select
      
      allocate( qin_chg_back(size( qin_chg )) )
      qin_chg_back(:) = 0._rp
      
      allocate( qlat_chg_back(size( qlat_chg(:,1)), size(qlat_chg(1,:)) ))
      qlat_chg_back(:,:) = 0._rp
      
!       if (allocated(bc%hyd_lat%x)) then
!          allocate(dof_back%qlat(1, size(bc%hyd_lat%x)))
!       else
!          allocate(dof_back%qlat(1, 1))
!       end if
!       if (allocated(bc%hyd_lat%t)) then
!          allocate(bc_back%hyd_lat%t(size(bc%hyd_lat%t)))
!       else
!          allocate(bc_back%hyd_lat%t(1))
!       end if
!       if (allocated(bc%hyd_lat%q)) then
!          allocate(bc_back%hyd_lat%q(size(bc%hyd_lat%q, 1), size(bc%hyd_lat%q, 2)))
!       else
!          allocate(bc_back%hyd_lat%q(1, 1))
!       end if
!       if (allocated(bc%hyd_lat%qsim)) then
!          allocate(bc_back%hyd_lat%qsim(size(bc%hyd_lat%qsim)))
!       else
!          allocate(bc_back%hyd_lat%qsim(1))
!       end if

         
#endif

   END SUBROUTINE alloc_back_vars


END MODULE m_tap_vars
