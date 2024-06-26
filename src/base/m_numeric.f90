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
!> \file m_numeric.f90
!! \brief This file includes m_numeric module.
!! \details The file includes only m_numeric module (see doc m_numeric module).



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module m_numeric
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> Module m_numeric.
!!
!! \details This module includes all defintion of numerical functions used in this softfaware.
MODULE m_numeric

   USE m_common
   USE m_linear_algebra
   USE m_mesh
!    USE  m_time_screen                                                                                              !NOADJ
   implicit none



CONTAINS



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Linear Interpolation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Linear interpolation
   !!
   !! \details From two real arrays  \f$tx\f$ and \f$ty\f$ such that \f$tx^i\f$-> \f$tx^i\f$ for \f$\forall i\f$, this function
   !! gives the value \f$y\f$ of a \f$x\f$ by linear interpolation.
   !! \param[in] tx Array \f$tx\f$.
   !! \param[in] ty Array \f$ty\f$.
   !! \param[in] x Value of  \f$x\f$.
   !! \return Value of \f$y\f$ 
   real(rp) FUNCTION linear_interp( tx , ty , x )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), dimension(:), intent(in)  ::  tx , ty
      real(rp)              , intent(in)  ::  x

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      real(rp)  ::  alpha1

      integer(ip)  ::  ind , imin, imax

      !================================================================================================================!
      !  Begin Function
      !================================================================================================================!

      imin = lbound( tx(:) , dim=1 )
      imax = ubound( tx(:) , dim=1 )
      
      if (imin == imax) then
      
         linear_interp = ty( imin )
      
      else if      ( x + zerom < tx( imin ) ) then

         alpha1 = ( x              - tx( imin     ) ) / &
                 ( tx( imin + 1 ) - tx( imin     ) )
                 
!          print *, "x + zerom < tx( imin ):", x, tx( imin     )

         linear_interp = ( one - alpha1 ) * ty( imin     ) + &
                                 alpha1   * ty( imin + 1 )

      else if ( x + zerom < tx( imax ) ) then

         ind = imin + 1
         
         do while ( tx( ind ) < x + zerom ) ; ind = ind + 1 ; end do

         alpha1 = ( x         - tx( ind - 1 ) ) / &
                 ( tx( ind ) - tx( ind - 1 ) )

         linear_interp = ( one - alpha1 ) * ty( ind - 1 ) + &
                                 alpha1   * ty( ind     )

      else

         alpha1 = ( x              - tx( imax - 1 ) ) / &
                 ( tx( imax     ) - tx( imax - 1 ) )

         linear_interp = ( one - alpha1 ) * ty( imax - 1 ) + &
                                 alpha1   * ty( imax     )

      end if

   END FUNCTION linear_interp




!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  P0 Interpolation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> constante interpolation
   !!
   !! \details From two real arrays  \f$tx\f$ and \f$ty\f$ such that \f$tx^i\f$-> \f$tx^i\f$ for \f$\forall i\f$, this function
   !! gives the value \f$y\f$ of a \f$x\f$ by linear interpolation.
   !! \param[in] tx Array \f$tx\f$.
   !! \param[in] ty Array \f$ty\f$.
   !! \param[in] x Value of  \f$x\f$.
   !! \return Value of \f$y\f$ 
   real(rp) FUNCTION p0_interp( tx , ty , x )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), dimension(:), intent(in)  ::  tx , ty
      real(rp)              , intent(in)  ::  x

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      integer(ip)  ::  ind , imin, imax

      !================================================================================================================!
      !  Begin Function
      !================================================================================================================!

      imin = lbound( tx(:) , dim=1 )
      imax = ubound( tx(:) , dim=1 )

      if      ( x + zerom < tx( imin ) ) then

         p0_interp = ty(imin)

      else if ( x + zerom < tx( imax ) ) then

         ind = imin + 1

         do while ( tx( ind ) < x + zerom ) ; ind = ind + 1 ; end do

         p0_interp = ty(ind)

         !if (x-tx(ind-1)<tx(ind)-x) then 

         !   p0_interp = ty(ind-1)

         !else 

         !  p0_interp = ty(ind)

         !end if

      else

         p0_interp = ty(imax)

      end if

   END FUNCTION p0_interp




!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
! Fourier Serie
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Computation of the Fourier serie.
   !!
   !! \details This function computes the Fourier series \f$Q_{in}(t)\f$ from Fourier coefficient and parameter \f$N_{FS}\f$.
   !! Fourier serie is given by:
   !!
   !! \f$Q_{in}(t)=\frac{a_{0}}{2}+\sum_{n=1}^{N_{FS}}\left(a_{n}\cos(nt\frac{2\pi}{T})+b_{n}\sin(nt\frac{2\pi}{T})\right)\f$ 
   !! where \f$\{a_{0};a_{n},b_{n}\} \f$ , \f$n\in[1..N_{FS}]\f$, are the Fourier coefficients and \f$T\f$ is the total 
   !! simulation time. The lower frequency represented by the Fourier serie is \f$1/T\f$ and the highest one 
   !! is \f$N_{FS}/T\f$. 
   !! \param[in] t Time of the simulation \f$t\f$.
   !! \param[in] a0 First indexe of the Fourier serie \f$a_{0}\f$.
   !! \param[in] n Number of frequency  \f$N_{FS}\f$.
   !! \param[in] A Fourier serie coefficient \f$a_{i}\f$.
   !! \param[in] B Fourier serie coefficient \f$b_{i}\f$.
   !! \return Value of \f$Q_{in}(t)\f$
   real(rp) FUNCTION Fourier_Serie( t,a0,n, A, B)

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!
      real(rp)                 , intent(in)  ::  t

      real(rp)                 , intent(in)  ::  a0

      real(rp)   , dimension(:), intent(in)  ::   A , B

      integer(ip), dimension(:), intent(in)  ::   n

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      integer(ip)  ::  ind , imin, imax,iloc

      !================================================================================================================!
      !  Begin Function
      !================================================================================================================!

      imin = lbound( A(:) , dim=1 )
      imax = ubound( A(:) , dim=1 )

      Fourier_Serie=a0*0.5_rp


      !do iloc=1,imax,1
      do iloc=1, size( A(:) )
         Fourier_Serie=Fourier_Serie+A(iloc)*cos((2._rp*n(iloc)*t*pi)/ts)+B(iloc)*sin((2._rp*n(iloc)*t*pi)/ts)
      end do


   END FUNCTION Fourier_Serie
                                                                                                               !<NOADJ


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  minmod operateur
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!
   !> Minmod operator.
   !!
   !! \details This function return the minmod between two real \f$x\f$ and \f$y\f$. 
   !! To recall:
   !!
   !!     if (x>=0 and y>=0):
   !!       return min(x,y)
   !!     if (x<=0 and y<=0):
   !!       return max(x,y)
   !!     else 
   !!       return 0
   !!
   !! \param[in] x Real \f$x\f$.
   !! \param[in] y Real \f$y\f$.
   !! \return mimmod computed.
   real(rp) FUNCTION minmod(x,y)
 
      implicit none

      real(rp)    :: x, y

      if ((x.ge.0).and.(y.ge.0)) then
         minmod=min(x,y)

      elseif ((x.le.0).and.(y.le.0)) then
         minmod=max(x,y)
      
      else 
         minmod=0

      endif
   END FUNCTION minmod


   !> Simpson integration rule.
   !!
   !! \details This function return integration by simpson method of a function \f$f\f$ between 
   !!\f$x_1\f$ and \f$x_2\f$.
   !! \param[in] f Function \f$f\f$ to integrate.
   !! \param[in] xmin Superior boundary \f$x_1\f$.
   !! \param[in] xmax Inferior boundary \f$x_2\f$.
   !! \param[in] nbp number of discretization step.
   !! \return mimmod computed.
   real(rp) FUNCTION Simpson( f , xmin , xmax , nbp )

      USE m_common

      implicit none

      real(rp)  ::  f  ;  external f

      real(rp), intent(in)  ::  xmin , xmax

      integer(ip), intent(in)  ::  nbp

      real(rp)  ::  step , sum_s , x_stage

      integer(ip)  ::  i_stage , nb_stage

      if ( modulo( nbp , 2 ) == 0 ) then

         nb_stage = nbp

      else

         nb_stage = nbp + 1

      end if

      step  =  ( xmax - xmin ) / real( nb_stage , 8 )

      sum_s  =  f( xmin )

      do i_stage = 1,nb_stage,2

         x_stage  =  xmin  +  real( i_stage , 8 ) * step

         sum_s  =  sum_s  +  4._rp * f( x_stage )

      end do

      do i_stage = 2,nb_stage-1,2

         x_stage  =  xmin  +  real( i_stage , 8 ) * step

         sum_s  =  sum_s  +  2._rp * f( x_stage )

      end do

      sum_s  =  sum_s  +  f( xmax )

      Simpson = step * sum_s / 3._rp

   END FUNCTION Simpson

   !> Find zero of a function.
   !!
   !! \details This searches for a zero of a function Func(X) between the given values B and C
   !! until the width of the interval (B,C) has collapsed to within a tolerance specified by
   !! the stopping criterion, ABS(B-C) <= 2.*(RW*ABS(B)+AE).The method used is an efficient
   !! combination of bisection and the secant rule.
   !! \author  L F Shampine and H A Watts
   !! \param[in] Func Function Name of the real valued external function.  This name
   !! must be in an EXTERNAL statement in the calling program.  F must be a function of one 
   !! real argument.
   !! \param[in] B One end of the interval (B,C).  The value returned for B usually is the better
   !! approximation to a zero of F.
   !! \param[in] C The other end of the interval (B,C).
   !! \param[in] RE Relative error used for RW in the stopping criterion. If the requested RE is 
   !! less than machine precision, then RW is set to approximately machine precision.
   !! \param[in] AE Absolute error used in the stopping criterion.  If the  given interval (B,C) contains
   !! the origin, then a nonzero value should be chosen for AE.
   !! \param[in] IFLAG A status code.  User must check IFLAG after each call. Control returns to the user
   !! from FZERO in all cases. XERROR does not process diagnostics in these cases.
   !!    - if ==1: B is within the requested tolerance of a zero. The interval (B,C) collapsed to the requested
   !! tolerance, the function changes sign in (B,C), and Func(X) decreased in magnitude as (B,C) collapsed.
   !!    - if ==2: Func(B) = 0. However, the interval (B,C) may not have collapsed to the requested tolerance.
   !!    - if ==3: B may be near a singular point of Func(X). The interval (B,C) collapsed to the requested 
   !! tolerance and the function changes sign in (B,C), but Func(X) increased in magnitude as (B,C) collapsed,i.e.
   !! abs(Func(B out)) > max(abs(Func(B in)),abs(Func(C in))).
   !!    - if ==4: No change in sign of Func(X) was found although the interval (B,C) collapsed to the requested 
   !! tolerance. The user must examine this case and decide whether  B is near a local minimum of Func(X), or 
   !! B is near a zero of even multiplicity, or neither of these.
   !!    - if ==5: Too many function evaluations used.


   RECURSIVE SUBROUTINE fzero( Func,b,c,re,ae,iflag )

   !***DATE WRITTEN   700901   (YYMMDD)
   !***REVISION DATE  861211   (YYMMDD)
   !***REVISION DATE  020410   (YYMMDD) [CONVERSION TO Fortran 90]
   !***CATEGORY NO.  F1B
   !***KEYWORDS  LIBRARY=SLATEC,TYPE=SINGLE PRECISION(FZERO-S DFZERO-D),
   !             BISECTION,NONLINEAR,NONLINEAR EQUATIONS,ROOTS,ZEROES,
   !             ZEROS
   !***AUTHOR  SHAMPINE,L.F.,SNLA
   !           WATTS,H.A.,SNLA
   !***PURPOSE  FZERO searches for a zero of a function Func(X) in a given
   !            interval (B,C).  It is designed primarily for problems
   !            where Func(B) and Func(C) have opposite signs.
   !***DESCRIPTION
   !
   !     Based on a method by T J Dekker
   !     written by L F Shampine and H A Watts
   !
   !            FZERO searches for a zero of a function Func(X) between
   !            the given values B and C until the width of the interval
   !            (B,C) has collapsed to within a tolerance specified by
   !            the stopping criterion, ABS(B-C) <= 2.*(RW*ABS(B)+AE).
   !            The method used is an efficient combination of bisection
   !            and the secant rule.
   !
   !     Description Of Arguments
   !
   !     Func,B,C,R,RE and AE are input parameters
   !     B,C and IFLAG are output parameters (flagged by an * below)
   !
   !        Func  - Name of the real valued external function.  This name
   !                must be in an EXTERNAL statement in the calling
   !                program.  F must be a function of one real argument.
   !
   !       *B     - One end of the interval (B,C).  The value returned for
   !                B usually is the better approximation to a zero of F.
   !
   !       *C     - The other end of the interval (B,C)
   !
   !        R     - A (better) guess of a zero of Func which could help in  **** Removed ***
   !                speeding up convergence.  If Func(B) and Func(R) have
   !                opposite signs, a root will be found in the interval
   !                (B,R); if not, but Func(R) and Func(C) have opposite
   !                signs, a root will be found in the interval (R,C);
   !                otherwise, the interval (B,C) will be searched for a
   !                possible root.  When no better guess is known, it is
   !                recommended that r be set to B or C; because if R is
   !                not interior to the interval (B,C), it will be ignored.
   !
   !        RE    - Relative error used for RW in the stopping criterion.
   !                If the requested RE is less than machine precision,
   !                then RW is set to approximately machine precision.
   !
   !        AE    - Absolute error used in the stopping criterion.  If the
   !                given interval (B,C) contains the origin, then a
   !                nonzero value should be chosen for AE.
   !
   !       *IFLAG - A status code.  User must check IFLAG after each call.
   !                Control returns to the user from FZERO in all cases.
   !                XERROR does not process diagnostics in these cases.
   !
   !                1  B is within the requested tolerance of a zero.
   !                   The interval (B,C) collapsed to the requested
   !                   tolerance, the function changes sign in (B,C), and
   !                   Func(X) decreased in magnitude as (B,C) collapsed.
   !
   !                2  Func(B) = 0.  However, the interval (B,C) may not have
   !                   collapsed to the requested tolerance.
   !
   !                3  B may be near a singular point of Func(X).
   !                   The interval (B,C) collapsed to the requested tol-
   !                   erance and the function changes sign in (B,C), but
   !                   Func(X) increased in magnitude as (B,C) collapsed,i.e.
   !                     abs(Func(B out)) > max(abs(Func(B in)),abs(Func(C in)))
   !
   !                4  No change in sign of Func(X) was found although the
   !                   interval (B,C) collapsed to the requested tolerance.
   !                   The user must examine this case and decide whether
   !                   B is near a local minimum of Func(X), or B is near a
   !                   zero of even multiplicity, or neither of these.
   !
   !                5  Too many function evaluations used.
   !
   !***REFERENCES  L. F. SHAMPINE AND H. A. WATTS, *FZERO, A ROOT-SOLVING
   !                 CODE*, SC-TM-70-631, SEPTEMBER 1970.
   !               T. J. DEKKER, *FINDING A ZERO BY MEANS OF SUCCESSIVE
   !                 LINEAR INTERPOLATION*, CONSTRUCTIVE ASPECTS OF THE
   !                 FUNDAMENTAL THEOREM OF ALGEBRA, EDITED BY B. DEJON
   !                 P. HENRICI, 1969.

   !------ Modifed by D. Henn
   !       1) Fortran 90 style
   !       2) Set immediate error if bracketing values Func(B) and Func(C) are not
   !          opposite signs

   IMPLICIT NONE

   REAL(rp),    INTENT( INOUT ) :: b, c
   REAL(rp),    INTENT( IN    ) :: re, ae
   INTEGER(ip), INTENT( OUT   ) :: iflag

   INTEGER(ip), PARAMETER :: MAX_ITER = 100

   REAL(rp)    a, acbs, acmb, aw, cmb, er, p, q
   REAL(rp)    fa, fb, fc, fx
   REAL(rp)    rw, tol
   INTEGER(ip) ic, count

   REAL(rp), EXTERNAL :: Func

   !------ Set error tolerances

   er = 2.0 * TINY(1.0)
   rw = MAX(re,er)
   aw = MAX(ae,0.0)

   count = 0 !Iteration count
   ic    = 0

   !------ Get function at bracketing values

   fc = Func( c )
   fb = Func( b )

   !------ Check if sign changes in initial interval [b,c]; if not set error

   IF( fc*fb > 0. )THEN
     iflag = 4
     GOTO 9999
   END IF

   a    = c
   fa   = fc
   acbs = ABS(b-c)
   fx   = MAX(ABS(fb),ABS(fc))

   !------ Iteration loop

   DO

   !------ Interchange such that |Func(b)| < |Func(c)|

     IF( ABS(fc) < ABS(fb) )THEN
       a  = b
       fa = fb
       b  = c
       fb = fc
       c  = a
       fc = fa
     END IF

     cmb  = 0.5*(c-b)
     acmb = ABS(cmb)
     tol  = rw*ABS(b) + aw

   !------ Test stopping criteria; set iflag appropriately

     IF( acmb <= tol )THEN
       IF( fb*fc > 0. )THEN
         iflag = 4
       ELSE IF( ABS(fb) > fx )THEN
         iflag = 3
       ELSE
         iflag = 1
       END IF
       IF( ABS(fb) > ABS(fc) )THEN  !Make sure last call is for smallest function value
         fc = Func( c )
       ELSE
         fb = Func( b )
       END IF
       EXIT

     ELSE IF( fb == 0.0 )THEN
      iflag = 2;  EXIT

     ELSE IF( count >= MAX_ITER )THEN
      iflag = 5; EXIT

     END IF

   !------ Calculate new iterate implicitly as b+p/q; arrange such that p > 0.

     p = (b-a)*fb
     q = fa-fb
     IF( p < 0. )THEN
       p = -p
       q = -q
     END IF

   !------ Update a and check for satisfactory reduction in the size of the
   !       bracketing interval every four increments; if not, perform bisection.

     a  = b
     fa = fb

     ic = ic + 1

     IF( ic >= 4 )THEN

       IF( 8.*acmb >= acbs )THEN
         p = HUGE(1.) !To force bisection
       ELSE
         ic = 0
         acbs = acmb
       END IF

     END IF

   !------ Get new iterate

     IF( p < ABS(q)*tol )THEN !If too small, increment by tolerance
       b = b + SIGN(tol,cmb)
     ELSE IF( p < cmb*q )THEN !Use secant rule if root is between b and (c+b)/2
       b = b + p/q
     ELSE                     !Othewise, use bisection
       b = 0.5*(c+b)
     END IF

   !------ Evaluate function at new iterate b

     fb = Func( b )

   !------ Decide whether next step is interpolation or extrapolation

     IF( fb*fc > 0. )THEN
       c  = a
       fc = fa
     END IF

     count = count + 1

   END DO

   9999 CONTINUE

   RETURN
   END



                                                                                 !>NOADJ

END MODULE m_numeric
