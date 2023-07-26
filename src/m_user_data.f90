MODULE m_user_data

   USE m_mesh
   USE m_common
   USE m_model

   implicit none


CONTAINS

   real(rp) FUNCTION s0_user( x , y )

      implicit none

      real(rp), intent(in)  ::  x , y

      print *, "The default m_user_data source file does not implement s0_user. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION s0_user


   real(rp) FUNCTION h0_user( x , y )

      implicit none

      real(rp), intent(in)  ::  x , y

      print *, "The default m_user_data source file does not implement h0_user. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION h0_user


   real(rp) FUNCTION q0_user( x , y )

      implicit none

      real(rp), intent(in)  ::  x , y

      print *, "The default m_user_data source file does not implement q0_user. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION q0_user


   real(rp) FUNCTION inflow_user( t , x , y )

      implicit none

      real(rp), intent(in)  ::  t, x , y

      print *, "The default m_user_data source file does not implement inflow_user. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION inflow_user


   real(rp) FUNCTION outflow_user( t , x , y )

      implicit none

      real(rp), intent(in)  ::  t, x , y

      print *, "The default m_user_data source file does not implement outflow_user. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION outflow_user


   real(rp) FUNCTION s_exact( x , y , t )

      implicit none

      real(rp), intent(in)  ::  x , y, t

      print *, "The default m_user_data source file does not implement s_exact. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION s_exact


   real(rp) FUNCTION q_exact( x , y , t )

      implicit none

      real(rp), intent(in)  ::  x , y, t

      print *, "The default m_user_data source file does not implement q_exact. "//& 
               "Please provide the m_user_data.f90 file in your case directory"
      stop

   END FUNCTION q_exact

END MODULE m_user_data
