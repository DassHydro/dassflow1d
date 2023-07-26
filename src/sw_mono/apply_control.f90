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
!> @file preissmann_time_step.f90
!> @brief This file includes the subroutine for solving the Shallow-Water Equations using Preissmann scheme with LPI 
!> and double sweep methods.


!> \brief Solve the Shallow-Water Equations using Preissmann scheme with LPI and double sweep methods.
subroutine apply_control(ctrl, mdl)
    use m_common
    use m_control, only: Control
    use m_mesh, only: strickler_type_powerlaw_h
    use m_sw_mono, only: Model
    use m_numeric
    use m_linear_algebra
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> @param dof Unknowns of the model
    type(Control), intent(inout) :: ctrl
    !> @param mdl Model
    type(Model), intent(inout) :: mdl
    !=LOCAL VARIABLES==================================================================================================!
    ! True if bathy slopes need to be recomputed
    logical :: recompute_bathy_slopes
    ! Index of first value in control vector
    integer(ip) :: first_value
    ! Iterator
    integer(ip) :: i
    ! Index of control item
    integer(ip) :: iitem
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of segment
    integer(ip) :: iseg
    ! Index of strickler parameter
    integer(ip) :: ipar
    ! Index of control value
    integer(ip) :: ivalue
    ! Index of first value in control vector
    integer(ip) :: last_value
    ! Index of boundary condition
    integer(ip) :: ibc
    ! Index of inflow condition
    integer(ip) :: iic
    ! Unobserved depth
    real(rp) :: h0
    ! Actual control values (in case of change of variable)
    real(rp), dimension(size(ctrl%x)) :: ctrl_values
    
    if (ctrl%nitems == 0) return
    
    if (ctrl%Bdemi%n > 0) then
    
        ctrl_values(:) = ctrl%x0(:)
        call csr_x_vec_add(ctrl%Bdemi, ctrl%x, ctrl_values(:))
        
    else
    
        ctrl_values(:) = ctrl%x(:)
    
    end if
    
    ! Loop on control items
    do iitem = 1, size(ctrl%items)
            
        first_value = ctrl%items(iitem)%offset + 1
        last_value = ctrl%items(iitem)%offset + ctrl%items(iitem)%nx
    
        if (ctrl%items(iitem)%id(1:2) == "BC") then
        
            read(ctrl%items(iitem)%id(3:5), '(I3)') ibc
            
            if (mdl%bc(ibc)%id == "discharge") then
            
                do i = first_value, last_value
                    mdl%bc(ibc)%ts%y(i-first_value+1) = max(ctrl_values(i), mdl%qeps)
                end do

            else if (mdl%bc(ibc)%id == "elevation") then
            
                do i = first_value, last_value
                    mdl%bc(ibc)%ts%y(i-first_value+1) = ctrl_values(i)
                end do
            
            end if
    
        else if (ctrl%items(iitem)%id(1:2) == "IC") then
        
            read(ctrl%items(iitem)%id(3:5), '(I3)') iic
            
            do i = first_value, last_value
                mdl%ic(iic)%ts%y(i-first_value+1) = max(ctrl_values(i), 0.0_rp)
            end do
    
        else if (ctrl%items(iitem)%id == "BATHY") then
        
            recompute_bathy_slopes = .true.
            
            ivalue = first_value
            do iseg = 1, mdl%msh%nseg

                if (allocated(mdl%msh%seg(iseg)%bathy_field%x)) then
                
                    ! Update y values of spatial field which is defined on current segment
                    do i = 1, size(mdl%msh%seg(iseg)%bathy_field%x)
                        mdl%msh%seg(iseg)%bathy_field%y(i) = ctrl_values(ivalue)
                        ivalue = ivalue + 1
                    end do
                    
                else 
                    ! Update y values of spatial field which is defined on current segment
                    do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
                        mdl%msh%cs(ics)%bathy = ctrl_values(ivalue)
                        if (mdl%msh%cs(ics)%bathy > mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps) then
                            mdl%msh%cs(ics)%bathy = mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps
                        end if
                        call compute_levels_cum(mdl%msh%cs(ics))
                        ivalue = ivalue + 1
                    end do
                end if
            end do

            ! DEPRECATED
!             if (allocated(mdl%msh%bathy_field%y)) then
!             
!                 ivalue = first_value
!                 do i = 1, size(mdl%msh%bathy_field%x)
!                     mdl%msh%bathy_field%y(i) = ctrl_values(ivalue)
!                     ivalue = ivalue + 1
!                 end do
!                 
!             else
!             
!                 ivalue = first_value
!                 do iseg = 1, mdl%msh%nseg
!                     do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
!                         mdl%msh%cs(ics)%bathy = ctrl_values(ivalue)
!                         if (mdl%msh%cs(ics)%bathy > mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps) then
!                             mdl%msh%cs(ics)%bathy = mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps
!                         end if
!                         call compute_levels_cum(mdl%msh%cs(ics))
!                         ivalue = ivalue + 1
!                     end do
!                 end do
!             
!             end if
    
        else if (ctrl%items(iitem)%id == "A0") then
            ivalue = first_value
            do iseg = 1, mdl%msh%nseg
                if (allocated(mdl%msh%seg(iseg)%bathy_field%y)) then
                    call abort_solver("Cannot compute A0 from control with bathymetry fields")
                end if
                do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
                    h0 = ctrl_values(ivalue) / mdl%msh%cs(ics)%level_widths(1)
                    mdl%msh%cs(ics)%bathy = mdl%msh%cs(ics)%level_heights(1) - h0
                    if (mdl%msh%cs(ics)%bathy > mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps) then
                        mdl%msh%cs(ics)%bathy = mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps
                    end if
                    call compute_levels_cum(mdl%msh%cs(ics))
                    ivalue = ivalue + 1
                end do
            end do
            
            ! DEPRECATED
!             if (allocated(mdl%msh%bathy_field%y)) then
!                 ! TODO
! !                 ivalue = first_value
! !                 do i = 1, size(mdl%msh%bathy_field%x)
! !                     mdl%msh%bathy_field%y(i) = ctrl_values(ivalue)
! !                     ivalue = ivalue + 1
! !                 end do
!                 
!             else
!             
!                 ivalue = first_value
!                 do iseg = 1, mdl%msh%nseg
!                     do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
!                         h0 = ctrl_values(ivalue) / mdl%msh%cs(ics)%level_widths(1)
!                         mdl%msh%cs(ics)%bathy = mdl%msh%cs(ics)%level_heights(1) - h0
!                         if (mdl%msh%cs(ics)%bathy > mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps) then
!                             mdl%msh%cs(ics)%bathy = mdl%msh%cs(ics)%level_heights(1) - 1.1 * mdl%heps
!                         end if
!                         call compute_levels_cum(mdl%msh%cs(ics))
!                         ivalue = ivalue + 1
!                     end do
!                 end do
!             
!             end if
    
        else if (ctrl%items(iitem)%id(1:4) == "KPAR") then
        
            read(ctrl%items(iitem)%id(5:6), '(I2)') ipar
      
            ivalue = first_value
            do iseg = 1, mdl%msh%nseg

                if (allocated(mdl%msh%seg(iseg)%strickler_fields) .and. &
                    size(mdl%msh%seg(iseg)%strickler_fields) > 0) then
                    do i = 1, size(mdl%msh%seg(iseg)%strickler_fields(1)%x)
                        if (mdl%msh%strickler_type_code /= strickler_type_powerlaw_h) then
                            mdl%msh%seg(iseg)%strickler_fields(ipar)%y(i) = max(1.0, ctrl_values(ivalue))
                        else
                            mdl%msh%seg(iseg)%strickler_fields(ipar)%y(i) = ctrl_values(ivalue)
                        end if
                        ivalue = ivalue + 1
                    end do

                else 
                    do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
                        if (mdl%msh%strickler_type_code /= strickler_type_powerlaw_h) then
                            mdl%msh%cs(ics)%strickler_params(ipar) = max(1.0, ctrl_values(ivalue))
                        else
                            mdl%msh%cs(ics)%strickler_params(ipar) = ctrl_values(ivalue)
                        end if
                        ivalue = ivalue + 1
                    end do

                end if 

            end do
    
#ifdef FLOODPLAIN_MODEL
        else if (ctrl%items(iitem)%id == "AFP") then
        
            ivalue = first_value
            do iseg = 1, mdl%msh%nseg
                do ics = mdl%msh%seg(iseg)%first_cs, mdl%msh%seg(iseg)%last_cs
                    mdl%msh%cs(ics)%alphaFP = max(1e-8, ctrl_values(ivalue))
                    ivalue = ivalue + 1
                end do
            end do
#endif
            
        end if
    
    end do
    
    if (recompute_bathy_slopes) then
        call bathy_slopes(mdl%msh)
    end if
   
end subroutine apply_control
