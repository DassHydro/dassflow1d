    !> Compute background covariance matrix
    !!
    !! \details Compute the blocks of the background error covariance matrix
    SUBROUTINE background_error_covariance_matrix( msh )

      use m_model
      use m_linear_algebra
      
      implicit none
      
      !=============================================================================================================!
      !  Interface Variables
      !=============================================================================================================!

      type(Mesh), intent(in   )  ::  msh
         
      !=============================================================================================================!
      !  Local Variables
      !=============================================================================================================!
         
      integer :: n  , m , l                                     ! Size of the correlation matrix
      real(rp), dimension(:), allocatable :: x           ! Abscissae for computing correlation matrix
      real(rp), dimension(:,:), allocatable :: C         ! Correlation matrix (column-major format)
      character(len=lchar)  ::  file_name! for B and Bdemi hydrograph lat

      if ( c_manning == 1 ) then
      select case( K_model )
      case ('K(h)')
         select case( K_covariance )
         case('constant_K_covariance')
            allocate(Bdemi_K_alpha(K_params%nb_diff_K_reachs))
            allocate(Bdemi_K_beta(K_params%nb_diff_K_reachs))
            do i=1,K_params%nb_diff_K_reachs
               Bdemi_K_alpha(i) = sigma_alpha 
               Bdemi_K_beta(i) = sigma_beta
            end do
         
        case('K_gaussian_covariance')
           !=========================================================================================!
           !  Compute correlation matrix
           !=========================================================================================!
            n = K_params%nb_diff_K_reachs
            allocate( Bdemi_K_alpha_m ( n , n ) )
            allocate( Bdemi_K_beta_m ( n , n ) )
            
            !!!!!!!!!!!!!!!!!!!
            ! ALPHA
            !!!!!!!!!!!!!!!!!!!
            allocate( x( n ) )
            allocate( C( n , n ) )

            x( : ) = K_params%cs_start(:)
            call correlation_matrix( n , x , mu_alpha , C ) ! Compute correlation matrix
            Bdemi_K_alpha_m( : , : ) = sigma_alpha**2 * C !Compute B matrix as sigma**2 * C
            deallocate( x )
            deallocate( C )
            !======================================================================================!
            !  Save B matrix
            !======================================================================================!
            open(50, file='min/B_alpha.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_K_alpha_m( : , i )  !Save B matrix
            end do
            close(50)
           
            call cholesky_inplace( n , Bdemi_K_alpha_m ) !Compute B^{1/2} by cholesky decomposition

            open(50, file='min/Bdemi_alpha.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_K_alpha_m( : , i ) ! Save B^{1/2} matrix
            end do
            close(50)

            !!!!!!!!!!!!!!!!!!!
            ! BETA
            !!!!!!!!!!!!!!!!!!!
            allocate( x( n ) )
            allocate( C( n , n ) )

            x( : ) = K_params%cs_start(:)
            call correlation_matrix( n , x , mu_beta , C ) ! Compute correlation matrix
            Bdemi_K_beta_m( : , : ) = sigma_beta**2 * C !Compute B matrix as sigma**2 * C
            deallocate( x )
            deallocate( C )
            !======================================================================================!
            !  Save B matrix
            !======================================================================================!
            open(50, file='min/B_beta.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_K_beta_m( : , i )  !Save B matrix
            end do
            close(50)
           
            call cholesky_inplace( n , Bdemi_K_beta_m ) !Compute B^{1/2} by cholesky decomposition

            open(50, file='min/Bdemi_beta.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_K_beta_m( : , i ) ! Save B^{1/2} matrix
            end do
            close(50)


         end select
!!!!!!!!!!!!!!!
 
      case default
         print*, 'K_model not taken into account'
      end select
      end if
         
      if ( c_hydrograph == 1 ) then
         
         !=========================================================================================!
         !  Compute correlation matrix
         !=========================================================================================!
            n = size( bc%hyd%q(:) )
            allocate( Bdemi_hydrograph( n , n ) )
            allocate( x( n ) )
            allocate( C( n , n ) )
            
            x( : ) = bc%hyd%t( : )
            call correlation_matrix( n , x , mu_hydrograph , C )
            
            !======================================================================================!
            !  Compute B matrix as sigma**2 * C
            !======================================================================================!
            Bdemi_hydrograph( : , : ) = sigma_hydrograph**2 * C
            deallocate( x )
            deallocate( C )
                  
            !======================================================================================!
            !  Save B matrix
            !======================================================================================!
            open(50, file='min/B_hydrograph.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_hydrograph( : , i )
            end do
            close(50)
            
            !======================================================================================!
            !  Compute B^{1/2} by cholesky decomposition
            !======================================================================================!
            
            call cholesky_inplace( n , Bdemi_hydrograph )
                  
            !======================================================================================!
            ! Save B^{1/2} matrix
            !======================================================================================!
            open(50, file='min/Bdemi_hydrograph.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_hydrograph( : , i )
            end do
            close(50)
            
         end if

         if ( lat_inflow == '1' .and. c_hydrograph == 1 ) then

         !=========================================================================================!
         !  Compute correlation matrix
         !=========================================================================================!
            n = size( bc%hyd_lat%q,2)  ! to do : ecrire ac virgule 1, 2 number of time steps for hydrographs
            m = size( bc%hyd_lat%q,1)  ! number of hydrographs
            allocate( Bdemi_hydrograph_lat( n , n , m ) )

            
            do l = 1, m

	      sigma_hydrographlat = bc%hyd_lat%params(l,1)
	      mu_hydrographlat = bc%hyd_lat%params(l,2)

              allocate( x( n ) )
              allocate( C( n , n ) )

              x( : ) = bc%hyd_lat%t( : )
              call correlation_matrix( n , x , mu_hydrographlat , C )

	      !======================================================================================!
	      !  Compute B matrix as sigma**2 * C
	      !======================================================================================!

	      Bdemi_hydrograph_lat( : , : , l ) = sigma_hydrographlat**2 * C

	      deallocate( x )
	      deallocate( C )
		    
	      !======================================================================================!
	      !  Save B matrix
	      !======================================================================================!write(file_name,'(A,I3.3)') 'min/hydrograph_lat_' , ite_min
	      write(file_name,'(A,I3.3,A)') 'min/B_hydrograph_lat' , l, '.dat'
	      open(50, file=file_name)

	      do i = 1, n
		  write(50, '(10000(1X,E15.7))') Bdemi_hydrograph_lat( : , i , l )
	      end do
	      close(50)

	      !======================================================================================!
	      !  Compute B^{1/2} by cholesky decomposition
	      !======================================================================================!
	      call cholesky_inplace( n , Bdemi_hydrograph_lat( : , : , l ) )
	      !======================================================================================!
	      ! Save B^{1/2} matrix
	      !======================================================================================!
	      write(file_name,'(A,I3.3,A)') 'min/Bdemi_hydrograph_lat' , l, '.dat'
	      open(50, file=file_name)

	      do i = 1, n
		  write(50, '(10000(1X,E15.7))') Bdemi_hydrograph_lat( : , i , l )
	      end do
	      close(50)
            enddo

         end if
         
         
         if ( c_bathy == 1 ) then
         
            !================================================================================!
            !  Compute correlation matrix
            !================================================================================!
            
            n = nb_bathy_control_pts
            allocate( Bdemi_bathy( n , n ) )
            allocate( x( n ) )
            allocate( C( n , n ) )
            do i = 1, n
               x(i) = msh%cs(bathy_indices(i))%x
            end do
            call correlation_matrix( n , x , mu_bathy , C )
            
            !================================================================================!
            !  Compute B matrix as sigma**2 * C
            !================================================================================!
            
            Bdemi_bathy( : , : ) = sigma_bathy**2 * C
            deallocate( x )
            deallocate( C )
            
            !================================================================================!
            ! Save B matrix
            !================================================================================!
            
            open(50, file='min/B_bathy.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_bathy( : , i )
            end do
            close(50)
      
            !================================================================================!
            !  Compute B^{1/2} by cholesky decomposition
            !================================================================================!
            call cholesky_inplace( n , Bdemi_bathy )
            
            !================================================================================!
            ! Save B^{1/2} matrix
            !================================================================================!
            
            open(50, file='min/Bdemi_bathy.dat')
            do i = 1, n
                write(50, '(10000(1X,E15.7))') Bdemi_bathy( : , i )
            end do
            close(50)
            
         else if ( c_bathy == 3 ) then
         
            select case( bathymetry )
         
               case('parametric')
         
                  !================================================================================!
                  !  Compute correlation matrix
                  !================================================================================!
                  
                  n = nb_bathy_control_pts
                  allocate( Bdemi_bathy( n , n ) )
                  allocate( x( n ) )
                  allocate( C( n , n ) )
                  x( : ) = msh%cs( bathy_first : bathy_last : bathy_step )%x
                  call correlation_matrix( n , x , mu_bathy , C )
                  
                  !================================================================================!
                  !  Compute B matrix as sigma**2 * C
                  !================================================================================!
                  
                  Bdemi_bathy( : , : ) = sigma_bathy**2 * C
                  deallocate( x )
                  deallocate( C )
                  
                  !================================================================================!
                  ! Save B matrix
                  !================================================================================!
                  
                  open(50, file='min/B_bathy.dat')
                  do i = 1, n
                     write(50, '(10000(1X,E15.7))') Bdemi_bathy( : , i )
                  end do
                  close(50)
            
                  !================================================================================!
                  !  Compute B^{1/2} by cholesky decomposition
                  !================================================================================!
                  call cholesky_inplace( n , Bdemi_bathy )
                  
                  !================================================================================!
                  ! Save B^{1/2} matrix
                  !================================================================================!
                  
                  open(50, file='min/Bdemi_bathy.dat')
                  do i = 1, n
                     write(50, '(10000(1X,E15.7))') Bdemi_bathy( : , i )
                  end do
                  close(50)
         
               case default               
                  write(*,*) "[ ERROR ] 'bathymetry' must be set to 'parametric' for change of "//&
                             "variable with c_bathy=3"
                  stop
         
            end select
            
         end if
        
       END SUBROUTINE background_error_covariance_matrix

          
      !> Compute correlation matrix
      !!
      !! \details .
      SUBROUTINE correlation_matrix( n , x , mu, C )

        USE m_common 
         implicit none
         
         !=========================================================================================!
         !  Arguments
         !=========================================================================================!
         
         integer(ip), intent(in) :: n                       !> Size of the correlation matrix
         real(rp), dimension(n), intent(in) :: x            !> Abscissae for computing correlation matrix
         real(rp), intent(in) :: mu                         !> Correlation length
         real(rp), dimension(n,n), intent(out) :: C         !> Correlation matrix (column-major format)
         


         if( mu < 0.0_rp ) then         
            write(*,*) "[ ERROR ] correlation length 'mu' must be positive or zero"
            stop
            
         else if( mu < 1e-12 ) then         
            C( : , : ) = 0.0_rp            
            do i = 1 , n            
               C( i , i ) = 1.0_rp               
            end do
            
         else         
            do i = 1 , n            
               do j = 1 , n
                  ! tmp : pour eviter d'ecrire des nombres inferieurs a 1e-90 dans les .dat         
                  if ( abs(x( j ) - x( i ))/mu > 200 ) then
                     C(i,j) = 0.0_rp
                  else
                     C( i , j ) = exp( -abs( x( j ) - x( i ) ) / mu )   ! exp kernel
                  end if
              end do               
            end do
         
         end if

       END SUBROUTINE correlation_matrix


       
    !> Compute inverse of background error covariance matrices
    !!
       !! \details Compute the inverse of the blocks of the background error covariance matrix
       SUBROUTINE calc_inv_cov_matrix()
         use m_model
         use m_linear_algebra
         
         implicit none

         ! TO BE COMPLETED
         
       END SUBROUTINE calc_inv_cov_matrix
