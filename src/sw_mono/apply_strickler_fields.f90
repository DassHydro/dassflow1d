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
!> @file apply_strickler_fields.f90
!> @brief TODO

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  TODO
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> @brief TODO
subroutine apply_strickler_fields(msh)
    use m_common
    use m_mesh, only: strickler_type_constant
    use m_mesh, only: strickler_type_powerlaw_h
    use m_mesh, only: Mesh
    use m_mesh, only: K_constant
    use m_mesh, only: K_powerlaw_h
    use m_numeric, only: linear_interp
    use m_numeric, only: p0_interp
    implicit none
    type(Mesh), intent(inout) :: msh
    
    ! Index
    integer(ip) :: i
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of segment
    integer(ip) :: iseg
    ! Index of field
    integer(ip) :: ifield
    ! Index
    integer(ip) :: ii
    
    do iseg = 1, msh%nseg
        
        if (.not. allocated(msh%seg(iseg)%strickler_fields)) return
        if (size(msh%seg(iseg)%strickler_fields) == 0) return

        do ics = msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
        
            do ifield = 1, size(msh%seg(iseg)%strickler_fields)
            
                if (msh%seg(iseg)%strickler_fields(ifield)%interp == "pt") then
                    i = ics - msh%seg(iseg)%first_cs + 1
                    msh%cs(ics)%strickler_params(ifield) = msh%seg(iseg)%strickler_fields(ifield)%y(i)
                else if (msh%seg(iseg)%strickler_fields(ifield)%interp == "segment") then
                    msh%cs(ics)%strickler_params(ifield) = msh%seg(iseg)%strickler_fields(ifield)%y(1)
                else if (msh%seg(iseg)%strickler_fields(ifield)%interp == "linear") then
                    msh%cs(ics)%strickler_params(ifield) = linear_interp(msh%seg(iseg)%strickler_fields(ifield)%x, &
                                                                         msh%seg(iseg)%strickler_fields(ifield)%y, &
                                                                         msh%cs(ics)%x)
                else if (msh%seg(iseg)%strickler_fields(ifield)%interp == "constant" .or. &
                         msh%seg(iseg)%strickler_fields(ifield)%interp == "p0") then
                    msh%cs(ics)%strickler_params(ifield) = p0_interp(msh%seg(iseg)%strickler_fields(ifield)%x, &
                                                                     msh%seg(iseg)%strickler_fields(ifield)%y, &
                                                                     msh%cs(ics)%x)
!                     if (ifield == 1) print *, ics, msh%cs(ics)%x, msh%cs(ics)%strickler_params(ifield)
!                     if (ifield == 1 .and. ics == msh%seg(iseg)%last_cs) then
!                         print *, msh%seg(iseg)%strickler_fields(ifield)%x
!                     end if
                end if
            end do
!             if (msh%seg(iseg)%strickler_fields(1)%interp == "segment") then
! 
! 	        msh%cs(ics)%strickler_params(1) = msh%seg(iseg)%strickler_fields(1)%y(1)
!                 if (msh%strickler_type_code == strickler_type_powerlaw_h) then
!                     msh%cs(ics)%strickler_params(2) = msh%seg(iseg)%strickler_fields(2)%y(1)
!                 else
!                     msh%cs(ics)%strickler_params(2) = 0.0
!                 end if
! 
!             else 
! 
!                 if (msh%seg(iseg)%strickler_fields(1)%interp == "linear") then
!                     msh%cs(ics)%strickler_params(1) = linear_interp(msh%seg(iseg)%strickler_fields(1)%x, &
!                                                                     msh%seg(iseg)%strickler_fields(1)%y, &
!                                                                     msh%cs(ics)%x)
!                     if (msh%strickler_type_code == strickler_type_powerlaw_h) then
!                         msh%cs(ics)%strickler_params(2) = linear_interp(msh%seg(iseg)%strickler_fields(2)%x, &
!                                                                         msh%seg(iseg)%strickler_fields(2)%y, &
!                                                                         msh%cs(ics)%x)
!                     else
!                         msh%cs(ics)%strickler_params(2) = 0.0
!                     end if
!                 end if
! 
!                 if (msh%seg(iseg)%strickler_fields(1)%interp == "p0") then
!                     ii = ubound( msh%seg(iseg)%strickler_fields(1)%x(:) , dim=1 )
!                     if (iseg == 3 .and. msh%cs(ics)%x>msh%seg(iseg)%strickler_fields(1)%x(ii)) then 
!                     !artificiellement sur la partie aval du negro, le strickler du negro amont de la confluence
!                         !ii = ubound(msh%seg(iseg-1)%strickler_fields(1)%y,dim=1)
!                         msh%cs(ics)%strickler_params(1) = msh%seg(iseg-1)%strickler_fields(1)%y(1)
!                         if (msh%strickler_type_code == strickler_type_powerlaw_h) then
!                             msh%cs(ics)%strickler_params(2) = msh%seg(iseg-1)%strickler_fields(2)%y(1)
!                         else
!                             msh%cs(ics)%strickler_params(2) = 0.0
!                         end if
!                     else 
!                         msh%cs(ics)%strickler_params(1) = p0_interp(msh%seg(iseg)%strickler_fields(1)%x, &
!                                                                     msh%seg(iseg)%strickler_fields(1)%y, &
!                                                                     msh%cs(ics)%x)
!                         !if (msh%cs(ics)%x==msh%cs(msh%seg(iseg)%last_cs)%x .and. iseg == 1 ) then 
!                             !print *, "test ", msh%cs(ics)%x, msh%cs(ics)%strickler_params(1)
!                             !print *, "test1 ", msh%seg(iseg)%strickler_fields(1)%x(1)
!                             !print *, "test2 ", msh%seg(iseg)%strickler_fields(1)%y(1)
!                         !end if
!                         if (msh%strickler_type_code == strickler_type_powerlaw_h) then
!                             msh%cs(ics)%strickler_params(2) = p0_interp(msh%seg(iseg)%strickler_fields(2)%x, &
!                                                                         msh%seg(iseg)%strickler_fields(2)%y, &
!                                                                         msh%cs(ics)%x)
!                         else
!                             msh%cs(ics)%strickler_params(2) = 0.0
!                         end if
!                     end if 
!                 end if
! 
!             end if 

        end do
            
    end do


    !if (.not. allocated(msh%strickler_fields)) return
 ! 
 !   if (msh%strickler_fields(1)%interp == "segment") then
 !   
 !       do iseg = 1, msh%nseg
 !       
 !           do ics = msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
 !           
 !               msh%cs(ics)%strickler_params(1) = msh%strickler_fields(1)%y(iseg)
 !               if (msh%strickler_type_code == strickler_type_powerlaw_h) then
 !                   msh%cs(ics)%strickler_params(2) = msh%strickler_fields(2)%y(iseg)
 !               else
 !                   msh%cs(ics)%strickler_params(2) = 0.0
 !               end if
 !               
 !           end do
 !           
 !       end do
 !       
 !   else 
 !   
 !        do iseg = 1, msh%nseg
 !       
 !           do ics = msh%seg(iseg)%first_cs, msh%seg(iseg)%last_cs
 !           
 !               if (msh%strickler_fields(1)%interp == "linear") then
 !                   msh%cs(ics)%strickler_params(1) = linear_interp(msh%strickler_fields(1)%x, &
 !                                                                   msh%strickler_fields(1)%y, &
 !                                                                   msh%cs(ics)%x)
 !                   if (msh%strickler_type_code == strickler_type_powerlaw_h) then
 !                       msh%cs(ics)%strickler_params(2) = linear_interp(msh%strickler_fields(2)%x, &
 !                                                                       msh%strickler_fields(2)%y, &
 !                                                                       msh%cs(ics)%x)
 !                   else
 !                       msh%cs(ics)%strickler_params(2) = 0.0
 !                   end if
 !               end if
 !       
 !           end do
!
!        end do
!        
!    end if
    
end subroutine
