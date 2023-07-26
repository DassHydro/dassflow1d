MODULE m_user_data

   USE m_mesh
   USE m_common
   USE m_model

   implicit none


CONTAINS


   !===================================================================================================================!
   !  Initial flow area
   !===================================================================================================================!

   real(rp) FUNCTION s0_user( x , y )

      implicit none

      real(rp), intent(in)  ::  x , y

      if (x.le.12.5) then
         s0_user  =  1._rp
      else 
         s0_user = 0._rp
      end if

   END FUNCTION s0_user


   !===================================================================================================================!
   !  Initial water depth
   !===================================================================================================================!

   real(rp) FUNCTION h0_user( x , y )

      implicit none

      real(rp), intent(in)  ::  x , y

      if (x.le.12.5) then
         h0_user  =  1._rp
      else 
         h0_user = 0._rp
      end if

   END FUNCTION h0_user


   !===================================================================================================================!
   !  Initial discharge
   !===================================================================================================================!

   real(rp) FUNCTION q0_user( x , y )

      implicit none

      real(rp), intent(in)  ::  x , y

      q0_user  = 0.0_rp

   END FUNCTION q0_user


   !===================================================================================================================!
   !  Inflow boundary condition
   !===================================================================================================================!

   real(rp) FUNCTION inflow_user( t , x , y )

      implicit none

      real(rp), intent(in)  ::  t , x , y

      inflow_user  = 1.75_rp 

   END FUNCTION inflow_user


   !===================================================================================================================!
   !  Outflow boundary condition
   !===================================================================================================================!

   real(rp) FUNCTION outflow_user( t , x , y )

      implicit none

      real(rp), intent(in)  ::  t , x , y

      outflow_user  =  1.0_rp

   END FUNCTION outflow_user


   !===================================================================================================================!
   !  Exact flow area
   !===================================================================================================================!

   real(rp) FUNCTION s_exact( x , y , t )

      implicit none

      real(rp), intent(in)  ::  x , y , t

      real(rp)  ::  x_a, x_b
      
      real(rp) :: Sl, x_0,L
      
      sl=1._rp
      x_0=12.5_rp
      L=25_rp
   
      x_a=x_0-t*sqrt(g*Sl)
      x_B=x_0+2_rp*t*sqrt(g*Sl)

      if (x.le.x_a) then
         s_exact=sl
      else if (x.ge.x_b) then 
         s_exact=0_rp
      else 
         s_exact=(4._rp/9._rp)*(sqrt(g*sl)- (x-x_0)/(2_rp*t))**2 * (1/g)
      end if

   END FUNCTION s_exact


   !===================================================================================================================!
   !  Exact water height
   !===================================================================================================================!

   real(rp) FUNCTION h_exact( x , y , t )

      implicit none

      real(rp), intent(in)  ::  x , y , t

      real(rp)  ::  x_a, x_b
      
      real(rp) :: Sl, x_0,L

      
      sl=1._rp
      x_0=12.5_rp
      L=25_rp
   
      x_a=x_0-t*sqrt(g*Sl)
      x_B=x_0+2_rp*t*sqrt(g*Sl)

      if (x.le.x_a) then
         h_exact=sl
      else if (x.ge.x_b) then 
         h_exact=0_rp
      else 
         h_exact=(4._rp/9._rp)*(sqrt(g*sl)- (x-x_0)/(2_rp*t))**2 * (1/g)
      end if

   END FUNCTION h_exact


   !===================================================================================================================!
   !  Exact discharge
   !===================================================================================================================!

   real(rp) FUNCTION q_exact( x , y , t )

      implicit none

      real(rp), intent(in)  ::  x , y , t

      real(rp)  ::  x_a, x_b
      
      real(rp) :: Sl, x_0,L


      sl=1._rp
      x_0=12.5_rp
      L=25_rp
   
      x_a=x_0-t*sqrt(g*Sl)
      x_B=x_0+2_rp*t*sqrt(g*Sl)

      if (x.le.x_a) then
         q_exact=0_rp
      else if (x.ge.x_b) then 
         q_exact=0_rp
      else 
         q_exact= (8._rp/(27._rp*g))*(((x-x_0)/t)+sqrt(g*Sl)) *(sqrt(g*sl)- (x-x_0)/(2_rp*t))**2 
      end if

   END FUNCTION q_exact

END MODULE m_user_data
