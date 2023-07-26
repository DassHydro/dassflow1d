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
!> \file constant_functions.f90
!! \brief This file includes m_regul_cost_coeff module and calc_one_minus_Fr2 subroutine.
!! \details The file only includes m_regul_cost_coeff module.
!! \todo Rename file ? Reorganize sources ?



!> Module m_regul_cost_coeff.
!!
!! \details This module includes the definition of the weighing coefficients for the regularization terms in the cost function.
!! These coeffs are defined in a specific module to prevent them to be differenciated by Tapenade.
MODULE m_regul_cost_coeff

   implicit none

   logical :: init_bathy_regul  !< Used to initialize the weighing coeff in the cost function for the bathy regularization


CONTAINS

   !> Calculation of the weighing coeff in the cost function for the bathy regularization
   !!
   !! \details The cost function is the sum of 3 terms: 
   !!    - The sum of the diff of the innonvation vector.
   !!    - The bathymetry regularization term.
  !!    - The hydograph regularization term.
  !!     J = J_obs + J_bathy_regul + J_hydro_regul
  !!     This subroutine computes the weighing term of J_bathy_regul
   !!
   !! \param[in] cost Value of the cost function.
  !! \param[in] cost_part Array containing the 2 parts of the cost : the observations and the regularization.
  !! \param[inout] alpha_cost_bathy Weighing coeff
   SUBROUTINE calc_bathy_regul_weighing_coeff( cost, cost_part, alpha_cost_bathy )

     USE m_numeric
      
      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), intent(inout)  ::  alpha_cost_bathy
      real(rp), intent(in)  ::  cost
      real(rp)  ::  cost_part(3)


            
      if ( init_bathy_regul ) then
         alpha_cost_bathy = 0.3
         print*, 'alpha_init', alpha_cost_bathy  !NOADJ
      else
         alpha_cost_bathy = 0.3*(cost - cost_part(1))/cost_part(2)
         print*, 'alpha', alpha_cost_bathy      !NOADJ
      end if

    END SUBROUTINE calc_bathy_regul_weighing_coeff

  END MODULE m_regul_cost_coeff




  !**********************************************************************************************************************!
  !**********************************************************************************************************************!
  !
  !  Compute 1-Fr^2
  !
  !**********************************************************************************************************************!
  !**********************************************************************************************************************!
  
  !>  Compute 1-Fr^2
  !!
  !! \details This subroutine compute the value \f$ 1-Fr^2 \f$, where \f$ Fr \f$ is the Froude number
  !! \param[in] q Discharge
  !! \param[in] w Top width
  !! \param[in] s Flow area
  !! \param[out] one_minus_Fr2 Value of \f$ 1-Fr^2 \f$
  SUBROUTINE calc_one_minus_Fr2(q, w, s, one_minus_Fr2) 
    USE m_numeric
   implicit none

   real(rp), intent(in) :: q, w, s
   real(rp), intent(out) :: one_minus_Fr2
   real(rp), parameter :: g = 9.81

   ! Result
#ifdef GEN_TAP
   one_minus_Fr2 = 0
#else
   one_minus_Fr2 = 1 - q**2 * w / (g * s**3) 
#endif

 END SUBROUTINE calc_one_minus_Fr2