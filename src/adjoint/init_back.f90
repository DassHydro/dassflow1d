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
subroutine init_back(cost_back, mdl, mdl_back, ctrl, ctrl_back, obs, obs_back)
    use m_common
    use m_control
    use m_control_back
    use m_sw_mono, only: Model
    use m_sw_mono_back, only: Model_back
    use m_numeric
    use m_obs
    use m_obs_back
    implicit none
    !=ARGUMENTS========================================================================================================!
    !> Adjoint of cost
    real(rp), intent(inout) :: cost_back
    !> @param mdl Model
    type(Model), intent(inout) :: mdl
    !> @param mdl_back Adjoint of model
    type(Model_back), intent(inout) :: mdl_back
    !> @param ctrl Control
    type(Control), intent(inout) :: ctrl
    !> @param ctrl_back Adjoint of control
    type(Control_back), intent(inout) :: ctrl_back
    !> @param obs Observations
    type(Observations), intent(inout) :: obs
    !> @param obs_back Adjoint of Observations
    type(Observations_back), intent(inout) :: obs_back
    !=LOCAL VARIABLES==================================================================================================!
    ! Iterator
    integer(ip) :: i
    ! Index of boundary condition
    integer(ip) :: ibc
    ! Index of cross-section
    integer(ip) :: ics
    ! Index of inflow condition
    integer(ip) :: iic
    ! Index of segment
    integer(ip) :: iseg
    
    ! Set adjoint of cost to 1.0
    cost_back = 1.0_rp
    
    ! (Re)allocate adjoint of model
    if (allocated(mdl_back%bc)) deallocate(mdl_back%bc)
    allocate(mdl_back%bc(size(mdl%bc)))
    do ibc = 1, size(mdl%bc)
        if (allocated(mdl%bc(ibc)%ts%y)) then
            allocate(mdl_back%bc(ibc)%ts%y(size(mdl%bc(ibc)%ts%y)))
        else
            allocate(mdl_back%bc(ibc)%ts%y(1))
        end if
        mdl_back%bc(ibc)%ts%y(:) = 0.0_rp
    end do
    
    if (allocated(mdl_back%ic)) deallocate(mdl_back%ic)
    if (allocated(mdl%ic) .and. mdl%nic > 0) then
        allocate(mdl_back%ic(size(mdl%ic)))
        do iic = 1, size(mdl%ic)
            if (allocated(mdl%ic(iic)%ts%y)) then
                allocate(mdl_back%ic(iic)%ts%y(size(mdl%ic(iic)%ts%y)))
            else
                allocate(mdl_back%ic(iic)%ts%y(1))
            end if
            mdl_back%ic(iic)%ts%y(:) = 0.0_rp
        end do
    else
        ! Allocate at least one dummy InflowCondition for the adjoint (size(mdl%ic, 1) is used)
        if (mdl%nic == 0) then
            if (.not. allocated(mdl%ic)) then
                allocate(mdl%ic(1))
            end if
        end if
        allocate(mdl_back%ic(1))
        allocate(mdl_back%ic(1)%ts%y(1))
    end if
    
    ! (Re)allocate adjoint of mesh
    !allocate(mdl_back%msh)
    if (allocated(mdl_back%msh%cs)) deallocate(mdl_back%msh%cs)
    allocate(mdl_back%msh%cs(size(mdl%msh%cs)))
    do ics = 1, size(mdl%msh%cs)
!         if (allocated(mdl%msh%cs(ics)%level_heights)) then
!             allocate(mdl_back%msh%cs(ics)%level_heights(size(mdl%msh%cs(ics)%level_heights)))
!         else
!             allocate(mdl_back%msh%cs(ics)%level_heights(1))
!         end if
!         mdl_back%msh%cs(ics)%level_heights(:) = 0.0_rp
!         if (allocated(mdl%msh%cs(ics)%level_heights)) then
!             allocate(mdl_back%msh%cs(ics)%level_widths(size(mdl%msh%cs(ics)%level_widths)))
!         else
!             allocate(mdl_back%msh%cs(ics)%level_widths(1))
!         end if
!         mdl_back%msh%cs(ics)%level_widths(:) = 0.0_rp
        if (allocated(mdl%msh%cs(ics)%strickler_params)) then
            allocate(mdl_back%msh%cs(ics)%strickler_params(size(mdl%msh%cs(ics)%strickler_params)))
        else
            allocate(mdl_back%msh%cs(ics)%strickler_params(1))
        end if
        mdl_back%msh%cs(ics)%strickler_params(:) = 0.0_rp
        if (allocated(mdl%msh%cs(ics)%area_cum)) then
            allocate(mdl_back%msh%cs(ics)%area_cum(0:size(mdl%msh%cs(ics)%area_cum)))
        else
            allocate(mdl_back%msh%cs(ics)%area_cum(0:1))
        end if
        mdl_back%msh%cs(ics)%area_cum(:) = 0.0_rp
        if (allocated(mdl%msh%cs(ics)%perim_cum)) then
            allocate(mdl_back%msh%cs(ics)%perim_cum(0:size(mdl%msh%cs(ics)%perim_cum)))
        else
            allocate(mdl_back%msh%cs(ics)%perim_cum(0:1))
        end if
        mdl_back%msh%cs(ics)%perim_cum(:) = 0.0_rp
    end do
    ! DEPRECATED
!     if (allocated(mdl_back%msh%bathy_field%y)) deallocate(mdl_back%msh%bathy_field%y)
!     if (allocated(mdl%msh%bathy_field%x)) then
!         allocate(mdl_back%msh%bathy_field%y(size(mdl%msh%bathy_field%x)))
!     else
!         allocate(mdl_back%msh%bathy_field%y(1))
!     end if
!     mdl_back%msh%bathy_field%y(:) = 0.0_rp
!     if (allocated(mdl_back%msh%strickler_fields)) then
!         do i = 1, size(mdl_back%msh%strickler_fields)
!             if (allocated(mdl_back%msh%strickler_fields(i)%y)) deallocate(mdl_back%msh%strickler_fields(i)%y)
!         end do
!         deallocate(mdl_back%msh%strickler_fields)
!     end if
!     if (allocated(mdl%msh%strickler_fields)) then
!         allocate(mdl_back%msh%strickler_fields(size(mdl%msh%strickler_fields)))
!         do i = 1, size(mdl_back%msh%strickler_fields)
!             allocate(mdl_back%msh%strickler_fields(i)%y(size(mdl%msh%strickler_fields(i)%y)))
!         end do
!     else
!         allocate(mdl_back%msh%strickler_fields(1))
!         allocate(mdl_back%msh%strickler_fields(1)%y(1))
!     end if
    if (allocated(mdl_back%msh%seg)) deallocate(mdl_back%msh%seg)
    allocate(mdl_back%msh%seg(size(mdl%msh%seg)))
    do iseg = 1, size(mdl%msh%seg)
        if (allocated(mdl_back%msh%seg(iseg)%bathy_field%y)) then
            if (allocated(mdl_back%msh%seg(iseg)%bathy_field%y)) deallocate(mdl_back%msh%seg(iseg)%bathy_field%y)
        end if
        if (allocated(mdl%msh%seg(iseg)%bathy_field%x)) then
            allocate(mdl_back%msh%seg(iseg)%bathy_field%y(size(mdl%msh%seg(iseg)%bathy_field%y)))
            mdl_back%msh%seg(iseg)%bathy_field%y(:) = 0.0_rp
        else
            allocate(mdl_back%msh%seg(iseg)%bathy_field%y(1))
            mdl_back%msh%seg(iseg)%bathy_field%y(1) = 0.0_rp
        end if
        if (allocated(mdl_back%msh%seg(iseg)%strickler_fields)) then
            do i = 1, size(mdl_back%msh%seg(iseg)%strickler_fields)
                if (allocated(mdl_back%msh%seg(iseg)%strickler_fields(i)%y)) deallocate(mdl_back%msh%seg(iseg)%strickler_fields(i)%y)
            end do
            deallocate(mdl_back%msh%seg(iseg)%strickler_fields)
        end if
        if (allocated(mdl%msh%seg(iseg)%strickler_fields)) then
            allocate(mdl_back%msh%seg(iseg)%strickler_fields(size(mdl%msh%seg(iseg)%strickler_fields)))
            do i = 1, size(mdl_back%msh%seg(iseg)%strickler_fields)
                allocate(mdl_back%msh%seg(iseg)%strickler_fields(i)%y(size(mdl%msh%seg(iseg)%strickler_fields(i)%y)))
                mdl_back%msh%seg(iseg)%strickler_fields(i)%y(:) = 0.0_rp
            end do
        else
            allocate(mdl%msh%seg(iseg)%strickler_fields(0))
            allocate(mdl_back%msh%seg(iseg)%strickler_fields(1))
            allocate(mdl_back%msh%seg(iseg)%strickler_fields(1)%y(1))
            mdl_back%msh%seg(iseg)%strickler_fields(1)%y(1) = 0.0_rp
        end if
    end do

    
    ! (Re)allocate adjoint of implicit matrix
    if (allocated(mdl%imp%ga)) then
      if (allocated(mdl_back%imp%ga)) deallocate(mdl_back%imp%ga)
      if (allocated(mdl_back%imp%gb)) deallocate(mdl_back%imp%gb)
      if (allocated(mdl_back%imp%gc)) deallocate(mdl_back%imp%gc)
      if (allocated(mdl_back%imp%gd)) deallocate(mdl_back%imp%gd)
      if (allocated(mdl_back%imp%ge)) deallocate(mdl_back%imp%ge)
      if (allocated(mdl_back%imp%gf)) deallocate(mdl_back%imp%gf)
      if (allocated(mdl_back%imp%cr)) deallocate(mdl_back%imp%cr)
      if (allocated(mdl_back%imp%cs)) deallocate(mdl_back%imp%cs)
      if (allocated(mdl_back%imp%ct)) deallocate(mdl_back%imp%ct)
      if (allocated(mdl_back%imp%ta1)) deallocate(mdl_back%imp%ta1)
      if (allocated(mdl_back%imp%ta2)) deallocate(mdl_back%imp%ta2)
      if (allocated(mdl_back%imp%ta3)) deallocate(mdl_back%imp%ta3)
      if (allocated(mdl_back%imp%ta4)) deallocate(mdl_back%imp%ta4)
      if (allocated(mdl_back%imp%tb1)) deallocate(mdl_back%imp%tb1)
      if (allocated(mdl_back%imp%tb2)) deallocate(mdl_back%imp%tb2)
      if (allocated(mdl_back%imp%anz)) deallocate(mdl_back%imp%anz)
      if (allocated(mdl_back%imp%rhs)) deallocate(mdl_back%imp%rhs)
      allocate(mdl_back%imp%ga(size(mdl%imp%ga)))
      allocate(mdl_back%imp%gb(size(mdl%imp%gb)))
      allocate(mdl_back%imp%gc(size(mdl%imp%gc)))
      allocate(mdl_back%imp%gd(size(mdl%imp%gd)))
      allocate(mdl_back%imp%ge(size(mdl%imp%ge)))
      allocate(mdl_back%imp%gf(size(mdl%imp%gf)))
      allocate(mdl_back%imp%cr(size(mdl%imp%cr)))
      allocate(mdl_back%imp%cs(size(mdl%imp%cs)))
      allocate(mdl_back%imp%ct(size(mdl%imp%ct)))
      allocate(mdl_back%imp%ta1(size(mdl%imp%ta1)))
      allocate(mdl_back%imp%ta2(size(mdl%imp%ta2)))
      allocate(mdl_back%imp%ta3(size(mdl%imp%ta3)))
      allocate(mdl_back%imp%ta4(size(mdl%imp%ta4)))
      allocate(mdl_back%imp%tb1(size(mdl%imp%tb1)))
      allocate(mdl_back%imp%tb2(size(mdl%imp%tb2)))
      allocate(mdl_back%imp%anz(size(mdl%imp%anz)))
      allocate(mdl_back%imp%rhs(size(mdl%imp%rhs)))
    else
      allocate(mdl%imp%ga(1))
      allocate(mdl%imp%gb(1))
      allocate(mdl%imp%gc(1))
      allocate(mdl%imp%gd(1))
      allocate(mdl%imp%ge(1))
      allocate(mdl%imp%gf(1))
      allocate(mdl%imp%cr(1))
      allocate(mdl%imp%cs(1))
      allocate(mdl%imp%ct(1))
      allocate(mdl%imp%ta1(1))
      allocate(mdl%imp%ta2(1))
      allocate(mdl%imp%ta3(1))
      allocate(mdl%imp%ta4(1))
      allocate(mdl%imp%tb1(1))
      allocate(mdl%imp%tb2(1))
      allocate(mdl%imp%anz(1))
      allocate(mdl%imp%rhs(1))
      allocate(mdl_back%imp%ga(1))
      allocate(mdl_back%imp%gb(1))
      allocate(mdl_back%imp%gc(1))
      allocate(mdl_back%imp%gd(1))
      allocate(mdl_back%imp%ge(1))
      allocate(mdl_back%imp%gf(1))
      allocate(mdl_back%imp%cr(1))
      allocate(mdl_back%imp%cs(1))
      allocate(mdl_back%imp%ct(1))
      allocate(mdl_back%imp%ta1(1))
      allocate(mdl_back%imp%ta2(1))
      allocate(mdl_back%imp%ta3(1))
      allocate(mdl_back%imp%ta4(1))
      allocate(mdl_back%imp%tb1(1))
      allocate(mdl_back%imp%tb2(1))
      allocate(mdl_back%imp%anz(1))
      allocate(mdl_back%imp%rhs(1))
    end if
    mdl_back%imp%ga(:) = 0.0_rp
    mdl_back%imp%gb(:) = 0.0_rp
    mdl_back%imp%gc(:) = 0.0_rp
    mdl_back%imp%gd(:) = 0.0_rp
    mdl_back%imp%ge(:) = 0.0_rp
    mdl_back%imp%gf(:) = 0.0_rp
    mdl_back%imp%cr(:) = 0.0_rp
    mdl_back%imp%cs(:) = 0.0_rp
    mdl_back%imp%ct(:) = 0.0_rp
    mdl_back%imp%ta1(:) = 0.0_rp
    mdl_back%imp%ta2(:) = 0.0_rp
    mdl_back%imp%ta3(:) = 0.0_rp
    mdl_back%imp%ta4(:) = 0.0_rp
    mdl_back%imp%tb1(:) = 0.0_rp
    mdl_back%imp%tb2(:) = 0.0_rp
    mdl_back%imp%anz(:) = 0.0_rp
    mdl_back%imp%rhs(:) = 0.0_rp
    
    ! (Re)allocate adjoint of unknowns
    if (allocated(mdl_back%dof%a)) deallocate(mdl_back%dof%a)
    if (allocated(mdl_back%dof%h)) deallocate(mdl_back%dof%h)
    if (allocated(mdl_back%dof%q)) deallocate(mdl_back%dof%q)
    allocate(mdl_back%dof%a(size(mdl%dof%a)))
    allocate(mdl_back%dof%h(size(mdl%dof%h)))
    allocate(mdl_back%dof%q(size(mdl%dof%q)))
    mdl_back%dof%a(:) = 0.0_rp
    mdl_back%dof%h(:) = 0.0_rp
    mdl_back%dof%q(:) = 0.0_rp
    
    ! (Re)allocate adjoint of control
    if (allocated(ctrl_back%x)) deallocate(ctrl_back%x)
    allocate(ctrl_back%x(size(ctrl%x)))
    ctrl_back%x(:) = 0.0_rp
    
    ! (Re)allocate adjoint of observations
    if (allocated(obs_back%est)) deallocate(obs_back%est)
    allocate(obs_back%est(size(obs%est, 1), size(obs%est, 2)))
    obs_back%est(:, :) = 0.0_rp

end subroutine init_back
