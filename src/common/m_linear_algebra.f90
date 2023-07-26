!======================================================================================================================!
!
!                    DassFlow1D Version 2.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures:
!    - 1D Shallow Module (1D Shallow Water Model, Finite Volume/Finite Difference Methods), i.e. the present code.
!    - 2D Shallow Module (2D Shallow Water Model, Finite Volume Method)
!    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current main developers or scientific contributers are:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               K. Larnier (C.S Communication and Systems & INSA Toulouse).
!               J. Monnier (INSA Toulouse & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA Toulouse & Mathematics Institute of Toulouse IMT).
!               P.-A. Garambois (INSA Strasbourg & ICUBE).
!               L. Pujol (CNES & INSA Strasbourg & ICUBE).
!  and former other developers (P. Brisset, R. Madec, M. Honnorat and J. Marin).
!
!  Scientific Contact : jerome.monnier@insa-toulouse.fr
!  Technical  Contact : kevin.larnier@c-s.fr
!                       frederic.couderc@math.univ-toulouse.fr
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
!> @file m_linear_algebra.f90
!> @brief This file includes the m_linear_algebra module.

!> @brief Module m_linear_algebra.
MODULE m_linear_algebra
    use m_common
    implicit none

    !> 2D Vector type
    !!
    !! \details This type is used for the definition of vector of size 2 ( Coordinates \f$(x,y)\f$ for example).
    type vec2d
        !> Coordinate \f$x\f$ of the vector.
        real(rp)  ::  x
        !> Coordinate \f$y\f$ of the vector.
        real(rp)  ::  y
    end type vec2d

    !> Dense Matrix type
    type Matrix
        !> Data array (COLUMN MAJOR : column index, row index)
        real(rp), dimension(:, :), allocatable  ::  m
    end type Matrix


    !> CSR Matrix type
    type MatrixCsr
        !> Size
        integer(ip) :: n = 0
        !> Number of non-zero entries
        integer(ip) :: nnz = 0
        !> Array of row indices for non-zero entries
        integer(ip), dimension(:), allocatable  ::  irow
        !> Data array (non zeros entries)
        integer(ip), dimension(:), allocatable  ::  icol
        !> Data array (non zeros entries)
        real(rp), dimension(:), allocatable  ::  anz
    end type MatrixCsr


    !> Block Matrix type
    type MatrixBlock
        !> Number of blocks
        integer(ip) :: n = 0
        !> Block dense matrices
        type(Matrix), dimension(:), allocatable  ::  blocks
    end type MatrixBlock
   


    contains
#ifndef CPP_ADJ
    
    
    subroutine matrixblock_initialise(mat, n, sizes)
        implicit none
        integer(ip), intent(in) :: n
        integer(ip), dimension(:), intent(in) :: sizes
        type(MatrixBlock), intent(out) :: mat
        
        integer :: i

        allocate(mat%blocks(n))
        do i = 1, n
            allocate(mat%blocks(i)%m(sizes(i), sizes(i)))
        end do
        
    end subroutine
    
    
    subroutine matrixblock_copy(src, dst)
        implicit none
        type(MatrixBlock), intent(in) :: src
        type(MatrixBlock), intent(out) :: dst
        
        integer :: i

        allocate(dst%blocks(size(src%blocks)))
        do i = 1, size(src%blocks)
            allocate(dst%blocks(i)%m(size(src%blocks(i)%m, 1), size(src%blocks(i)%m, 2)))
            dst%blocks(i)%m(:, :) = src%blocks(i)%m(:, :)
        end do
        
    end subroutine
    
    
    subroutine matrixcsr_initialise(mat, n, nnz)
        implicit none
        integer(ip), intent(in) :: n
        integer(ip), intent(in) :: nnz
        type(MatrixCsr), intent(out) :: mat

        allocate(mat%irow(nnz))
        allocate(mat%icol(nnz))
        allocate(mat%anz(nnz))
        mat%irow(:) = 0
        mat%icol(:) = 0
        mat%anz(:) = 0.0_rp
        
    end subroutine
    
    
    subroutine matrixcsr_from_matrixblock(matblock, matcsr, threshold)
        implicit none
        type(MatrixBlock), intent(in) :: matblock
        real(rp), intent(in), optional :: threshold
        type(MatrixCsr), intent(inout) :: matcsr
        
        integer(ip) :: i
        integer(ip) :: icol
        integer(ip) :: inz
        integer(ip) :: irow
        integer(ip) :: offset
        real(rp) :: nz_threshold
        
        ! CHECK-UP
        do i = 1, size(matblock%blocks)
            if (size(matblock%blocks(i)%m, 1) /= size(matblock%blocks(i)%m, 2)) then
                call abort_solver("Blocks must be square matrices")
            end if
        end do
        
        ! Setup non-zero threshold
        if (present(threshold)) then
            nz_threshold = threshold
        else
            nz_threshold = 1e-15
        end if
        
        ! Count non-zero entries and total size
        matcsr%nnz = 0
        matcsr%n = 0
        do i = 1, size(matblock%blocks)
            matcsr%n = matcsr%n + size(matblock%blocks(i)%m, 1)
            do irow = 1, size(matblock%blocks(i)%m, 1)
                do icol = 1, size(matblock%blocks(i)%m, 2)
                    if (matblock%blocks(i)%m(irow, icol) >= nz_threshold) matcsr%nnz = matcsr%nnz + 1
                end do
            end do
        end do
        
        ! Set matrix arrays
        if (allocated(matcsr%irow)) deallocate(matcsr%irow)
        if (allocated(matcsr%icol)) deallocate(matcsr%icol)
        if (allocated(matcsr%anz)) deallocate(matcsr%anz)
        allocate(matcsr%irow(matcsr%nnz))
        allocate(matcsr%icol(matcsr%nnz))
        allocate(matcsr%anz(matcsr%nnz))
        inz = 0
        offset = 0
        do i = 1, size(matblock%blocks)
            do irow = 1, size(matblock%blocks(i)%m, 1)
                do icol = 1, size(matblock%blocks(i)%m, 2)
                    if (matblock%blocks(i)%m(irow, icol) >= nz_threshold) then
                        inz = inz + 1
                        matcsr%irow(inz) = offset + irow-1
                        matcsr%icol(inz) = offset + icol-1
                        matcsr%anz(inz) = matblock%blocks(i)%m(irow, icol)
                    end if
                end do
            end do
            offset = offset + size(matblock%blocks(i)%m, 1)
        end do
        
    end subroutine
    
    
    subroutine matrixcsr_from_numpy_array(array, mat, threshold)
        implicit none
        real(rp), dimension(:, :), intent(in) :: array
        real(rp), intent(in), optional :: threshold
        type(MatrixCsr), intent(out) :: mat
        
        integer(ip) :: icol
        integer(ip) :: inz
        integer(ip) :: irow
        real(rp) :: nz_threshold
        
        ! CHECK-UP
        if (size(array, 1) /= size(array, 2)) then
            call abort_solver("Array must be square")
        end if
        
        ! Setup non-zero threshold
        if (present(threshold)) then
            nz_threshold = threshold
        else
            nz_threshold = 1e-15
        end if
        
        ! Count non-zero entries
        mat%nnz = 0
        do irow = 1, size(array, 1)
            do icol = 1, size(array, 2)
                if (array(irow, icol) >= nz_threshold) mat%nnz = mat%nnz + 1
            end do
        end do
        
        ! Set matrix arrays
        mat%n = size(array, 1)
        allocate(mat%irow(mat%nnz))
        allocate(mat%icol(mat%nnz))
        allocate(mat%anz(mat%nnz))
        inz = 0
        do irow = 1, size(array, 1)
            do icol = 1, size(array, 2)
                if (array(irow, icol) >= nz_threshold) then
                    inz = inz + 1
                    mat%irow(inz) = irow-1
                    mat%icol(inz) = icol-1
                    mat%anz(inz) = array(irow, icol)
                end if
            end do
        end do
        
    end subroutine
    
! ! ! 
! ! ! 
! ! ! 
! ! !    !> Compute the Cholesky decomposition in place (lower triangular part is L, upper triangular part is L.T)
! ! !     subroutine cholesky_inplace_csr( A )
! ! !         implicit none
! ! !         !=ARGUMENTS====================================================================================================!
! ! !         !> Matrix A in CSR format
! ! !         type(MatrixCsr), intent(inout) :: A
! ! !         !=LOCAL VARIABLES==============================================================================================!
! ! !         integer(ip)  ::  i, j, k
! ! !         real(rp) :: sum1
! ! !         real(rp), dimension(n)  ::  diag
! ! !         
! ! !         inz = 1
! ! !         do i = 1 , n
! ! !           if (A%irow(inz) > i) cycle
! ! !           do while(A%irow(inz) == i)
! ! !           
! ! !               sum1 = A%anz(inz)
! ! !               
! ! !               ! TODO
! ! !               
! ! !               do k = i-1, 1 , -1
! ! !               
! ! !                 sum1 = sum1 - A( k , i ) * A( k , j )
! ! !                 
! ! !               end do
! ! !               
! ! !               if( i == j ) then
! ! !               
! ! !                 if( sum1 <= 0.0 ) then
! ! !                 
! ! !                     write(*,*) "[ ERROR ] Cholesky decomposition failed on sum1=", sum1
! ! !                     stop
! ! !                     
! ! !                 end if
! ! !                 
! ! !                 diag( i ) = sqrt( sum1 )
! ! !                 
! ! !               else
! ! !               
! ! !                 A( i , j ) = sum1 / diag( i )
! ! ! 
! ! !               end if
! ! ! 
! ! !           end do
! ! !           
! ! !         end do
! ! ! 
! ! !       !================================================================================================================!
! ! !       !  Apply computed diagonal values
! ! !       !================================================================================================================!
! ! !       
! ! !       do i = 1 , n
! ! !       
! ! !          A( i , i) = diag( i )
! ! !          
! ! !       end do
! ! ! ! 
! ! ! !       !================================================================================================================!
! ! ! !       !  Compute upper part by symmetry
! ! ! !       !================================================================================================================!
! ! ! !       
! ! ! !       do i = 1 , n
! ! ! !       
! ! ! !          do j = i+1 , n
! ! ! !       
! ! ! !             A( j , i ) = A( i , j )
! ! ! ! 
! ! ! !          end do
! ! ! !          
! ! ! !       end do
! ! ! 
! ! !       !================================================================================================================!
! ! !       !  Put zeros on upper part
! ! !       !================================================================================================================!
! ! !       
! ! !       do i = 1 , n
! ! !       
! ! !          do j = i+1 , n
! ! !       
! ! !             A( j , i ) = 0.0_rp
! ! ! 
! ! !          end do
! ! !          
! ! !       end do
! ! ! 
! ! !    end SUBROUTINE cholesky_inplace
    


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compute Cholesky decomposition for a dense matrix
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!
   !> Compute the Cholesky decomposition in place (lower triangular part is L, upper triangular part is L.T)
   !!
   !! \param[in] n Size of the Matrix A
   !! \param[inout] A Matrix A in dense format
   SUBROUTINE cholesky_inplace( n, A )

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      integer(ip), intent(in)  ::  n

      real(rp), dimension(n,n)  ::  A

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      integer(ip)  ::  i, j, k
      
      real(rp) :: sum1
      
      real(rp), dimension(n)  ::  diag

      !================================================================================================================!
      !  Compute diagonal values and lower trianguar part
      !================================================================================================================!

      do i = 1 , n
      
         do j = i , n
         
            sum1 = A( j , i )
            
            do k = i-1, 1 , -1
            
               sum1 = sum1 - A( k , i ) * A( k , j )
               
            end do
            
            if( i == j ) then
            
               if( sum1 <= 0.0 ) then
               
                  write(*,*) "[ ERROR ] Cholesky decomposition failed on sum1=", sum1
                  stop
                  
               end if
               
               diag( i ) = sqrt( sum1 )
               
            else
            
               A( i , j ) = sum1 / diag( i )

            end if

         end do
         
      end do

      !================================================================================================================!
      !  Apply computed diagonal values
      !================================================================================================================!
      
      do i = 1 , n
      
         A( i , i) = diag( i )
         
      end do
! 
!       !================================================================================================================!
!       !  Compute upper part by symmetry
!       !================================================================================================================!
!       
!       do i = 1 , n
!       
!          do j = i+1 , n
!       
!             A( j , i ) = A( i , j )
! 
!          end do
!          
!       end do

      !================================================================================================================!
      !  Put zeros on upper part
      !================================================================================================================!
      
      do i = 1 , n
      
         do j = i+1 , n
      
            A( j , i ) = 0.0_rp

         end do
         
      end do

   end SUBROUTINE cholesky_inplace
! 
! 
! !**********************************************************************************************************************!
! !**********************************************************************************************************************!
! !
! !  Compute Cholesky decomposition for a dense matrix
! !
! !**********************************************************************************************************************!
! !**********************************************************************************************************************!
!    !> Compute the Cholesky decomposition in place (lower triangular part is L, upper triangular part is L.T)
!    !!
!    !! \param[in] n Size of the Matrix A
!    !! \param[in] mat Matrix A in dense format
!    SUBROUTINE cholesky_inplace( n, mat )
! 
!       !================================================================================================================!
!       !  Interface Variables
!       !================================================================================================================!
! 
!       integer(ip), intent(in)  ::  n
! 
!       type(matrix)  ::  mat
! 
!       !================================================================================================================!
!       !  Interface Variables
!       !================================================================================================================!
! 
!       integer(ip)  ::  i, j
! 
!       !================================================================================================================!
!       !  Compute lower trianguar part
!       !================================================================================================================!
! 
!       mat%m(1,1) = sqrt(mat%m(1,1))
!       
!       do j = 2, n
!       
!          mat%m(1,j) = mat%m(j,1) / mat%m(1,1)
!          
!       end do
!       
!       do i = 2, n
!       
!          do k = 1, i-1
!          
!             mat%m(i,i) = mat%m(i,i) - mat%m(k,i)**2
!             
!          end do
!          
!          mat%m(i,i) = sqrt(mat%m(i,i))
!          
!          do j = i+1, n
!             
!             mat%m(i,j) = mat%m(j,i)
!             
!             do k = 1, i-1
!             
!                mat%m(i,j) = mat%m(i,j) - mat%m(k,i) * mat%m(k,j)
!                
!             end do
!             
!             mat%m(i,j) = mat%m(i,j) / mat%m(i,i)
!             
!          end do
! 
!       end do
! 
!       !================================================================================================================!
!       !  Compute upper triangular part (symmetricity)
!       !================================================================================================================!
! 
!       do i = 1, n
!       
!          do j = i+1, n
!       
!             mat%m(j,i) = mat%m(i,j)
!             
!          end do
!          
!       end do
! 
!    end SUBROUTINE cholesky_inplace


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Solve Ax=b where A is a Cholesky decomposition
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!
   !> Solve Ax=b using the Cholesky decomposition of A (lower triangular part is L, upper triangular part is L.T)
   !!
   !! \param[in] n Size of the Matrix A
   !! \param[in] mat Matrix A in dense format
   !! \param[in] b Right hand side vector
   !! \param[in] x Solution
   SUBROUTINE solve_using_cholesky_inplace( n, mat, b, x )

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      integer(ip), intent(in)  ::  n

      type(matrix), intent(in) ::  mat
      
      real(rp), dimension(n), intent(in)  ::  b
      
      real(rp), dimension(n), intent(out)  ::  x

      !================================================================================================================!
      !  
      !================================================================================================================!

      integer(ip)  ::  i, j
      real(rp), dimension(n)  ::  y

      !================================================================================================================!
      !  Step 1 : solve Ly = b
      !================================================================================================================!
      y(1) = b(1) / mat%m(1,1)
      do i = 2, n
      
         y(i) = b(i)
         
         do j = 1, i-1
         
            y(i) = y(i) - mat%m(j,i) * b(j)
         
         end do
         y(i) = y(i) / mat%m(i,i)
        
      end do

      !================================================================================================================!
      !  Step 2 : solve L'x = y
      !================================================================================================================!
      x(n) = y(n) / mat%m(n,n)
      do i = n-1, 1, -1
      
         x(i) = y(i)
         
         do j = i+1, n
         
            x(i) = x(i) - mat%m(j,i) * y(j)
         
         end do
         x(i) = x(i) / mat%m(i,i)
        
      end do

   end SUBROUTINE solve_using_cholesky_inplace
   
#endif

    !> Compute linear regression
    subroutine linear_regression(X, Y, a, b)
        implicit none
        !=ARGUMENTS====================================================================================================!
        !> Vector X
        real(rp), dimension(:), intent(in) :: X
        !> Vector Y
        real(rp), dimension(:), intent(in) :: Y
        !> Slope
        real(rp), intent(out) :: a
        !> Intercept
        real(rp), intent(out) :: b
        !=LOCAL VARIABLES==============================================================================================!
        ! Multiplier
        real(rp) :: mult
        ! Sum of values in X
        real(rp) :: sum_x
        ! Sum of values in Y
        real(rp) :: sum_y
        ! Sum of products xi * xi
        real(rp) :: sum_xx
        ! Sum of products xi * yi
        real(rp) :: sum_xy
        
        sum_x = sum(x)
        sum_y = sum(y)
        sum_xx = sum(x*x)
        sum_xy = sum(x*y)
        mult = 1.0 / (size(x) * sum_xx - sum_x * sum_x)
        a = (size(x) * sum_xy - sum_x * sum_y) * mult
        b = (sum_xx * sum_y - sum_xy * sum_x) * mult
        
    end subroutine

    !> Compute Y += A.X
    subroutine csr_x_vec_add(A, X, Y)
        implicit none
        !=ARGUMENTS====================================================================================================!
        !> Matrix A in CSR format
        type(MatrixCsr), intent(inout) :: A
        !> Vector X
        real(rp), dimension(:), intent(in) :: X
        !> Vector Y
        real(rp), dimension(:), intent(inout) :: Y
        !=LOCAL VARIABLES==============================================================================================!
        ! Index of column
        integer(ip) :: icol
        ! Index of non-zero entry
        integer(ip) :: inz
        ! Index of row
        integer(ip) :: irow
        
        do inz = 1, A%nnz
            irow = A%irow(inz)+1
            icol = A%icol(inz)+1
            Y(irow) = Y(irow) + A%anz(inz) * X(icol)
        end do
        
        
    end subroutine

end module m_linear_algebra
