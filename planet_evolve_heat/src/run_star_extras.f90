! ***********************************************************************
!
!   Copyright (C) 2011  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none
      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_energy => my_other_energy
      end subroutine extras_controls
      

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 12
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         use const_def, only: pi
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i, j, k, l, nz
         real(dp), allocatable :: q_c(:), scale_height(:), scale_height_t(:), dr(:)
         real(dp) :: irradiation_dq, xq
         real(dp) :: rho_mean, q_0, volume, b_scaling, F_factor

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         nz = s% nz
         allocate(q_c(nz),scale_height(nz),scale_height_t(nz),dr(nz))
         do l = 2, nz-1
            dr(l) = 0.5*(s%r(l-1) - s%r(l+1))
         end do
         dr(1) = dr(2)  
         dr(nz) = dr(nz-1) ! The last and first dr do not matter (we will not integrate over themS)

         !write(*, *) 'dr', dr

         irradiation_dq = pi4*s% r(1)*s% r(1)*s% column_depth_for_irradiation/s% xmstar
         do i = 1, s% nz, 1
            xq = 1 - s% q(i)
            if (irradiation_dq < xq) exit
         end do

         do j = 1, nz, 1
            if (s% peos(j) > 1e12) exit
         end do

         do k = 1, nz, 1
            if ((1 - s% q(k)) > 2*(1 - s% q(i))) exit
         end do

         names(1) = "conv_luminosity"
         vals(1) = s% l(k) / lsun ! in luminosities

         names(2) = "conv_temperature"
         vals(2) = s% t(i) ! in K

         names(3) = "conv_radius"
         vals(3) = s% r(i) ! in cm

         names(4) = "conv_mass"
         vals(4) = s% m(i) ! in g

         names(5) = "conv_density"
         vals(5) = s% rho(i) ! in g cm^-3

         names(6) = "dyn_luminosity"
         vals(6) = s% l(j) / lsun ! in luminosities

         names(7) = "dyn_temperature"
         vals(7) = s% t(j) ! in K

         names(8) = "dyn_radius"
         vals(8) = s% r(j) ! in cm

         names(9) = "dyn_mass"
         vals(9) = s% m(j) ! in g

         names(10) = "dyn_density"
         vals(10) = s% rho(j) ! in g cm^-3

         q_c = 2 * s%cp * s%t * s%rho**2 * s%conv_vel**3 / s%peos / s%chiT * s%chiRho ! convective flux
         q_0 = q_c(j) 
         ! q_0 = q_c(nz-1) * s%r(nz-1) * s%r(nz-1) / s%r(j) / s%r(j) ! reference convective flux at the outer boundary of the dyanmo
         ! write(*, *) 'q_0', q_0
         ! write(*, *) 'l/4 pi r^2', s%l(j) / (4 * pi * s%r(j)**2)
         volume = 4 * pi  * (s%r(j)**3 - s%r(nz-1)**3) / 3 ! dynamo shell volume
         rho_mean = (s%m(j) -  s%m(nz-1)) / volume

         do l = 1, nz, 1
            scale_height(l) = min(s%scale_height(l), s%r(j) - s%r(nz-1))
         end do
         scale_height_t = s%peos / (s%rho * s%grav * s%grada)
         !write(*, *) 'scale_height', scale_height(1:5)
         !write(*, *) 'scale_height_t', scale_height_t(1:5)

         !F_factor = sum((scale_height(j:nz-1)/scale_height_t(j:nz-1)*q_c(j:nz-1)/q_0)**(2d0/3d0) * (s%rho(j:nz-1)/rho_mean)**(1d0/3d0) * s%r(j:nz-1)**2 * dr(j:nz-1))
         !F_factor = 4 * pi * F_factor / volume ! This is actually F**(2/3)
         !write(*, *) 'F_factor_1', F_factor**(3d0/2d0)

         F_factor = sum((scale_height(j:nz-1)/scale_height_t(j:nz-1)*q_c(j:nz-1)/q_0)**(2d0/3d0) * (s%rho(j:nz-1)/rho_mean)**(1d0/3d0) * s%dm(j:nz-1) / s%rho(j:nz-1)) !s%dq(j:nz-1) * (s%m(j) -  s%m(nz-1))
         F_factor = F_factor / volume ! This is actually F**(2/3)
         write(*, *) 'F_factor:', F_factor**(3d0/2d0)

         names(11) = "F_factor"
         vals(11) = F_factor**(3d0/2d0)

         ! Christensen2009 scaling law:
         ! <B>**2/(2 mu0) = c fohm <rho>**(1/3) (F q0)**(2/3)
         ! where we set c = 0.63, fohm = 1
         b_scaling = 2 * 1.256637e-6 * 0.63 * 1 * (1e3*rho_mean)**(1.0/3.0) * F_factor * (q_0/1e3)**(2.0/3.0)
         b_scaling = sqrt(b_scaling)*1e4

         ! write(*, *) '2 mu0 c', 2 * 1.256637e-6 * 0.63
         ! write(*, *) 'F_factor', F_factor
         ! write(*, *) 'q0', q_0/1e3
         ! write(*, *) 'rho', 1e3*rho_mean
         ! write(*, *) 'b_scaling', b_scaling

         names(12) = "b_scaling"
         vals(12) = b_scaling

      end subroutine data_for_extra_history_columns
      
      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         real(dp) :: dgamma1(nz),b_crit(nz)
         real(dp) :: delta_grad,Q
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         do k = 2, nz-1
            dgamma1(k) = 2d0 * (s%gamma1(k+1) - s%gamma1(k)) / (s%gamma1(k+1) + s%gamma1(k))
            dgamma1(k) = dgamma1(k) * 2d0 * (s%peos(k) - s%peos(k-1)) / (s%peos(k) + s%peos(k-1))
         end do
         dgamma1(1)=0
         dgamma1(nz)=0

         do k = 1, nz    
            delta_grad = s%gradr(k) - s%grada(k)
            Q = 1 + 4d0 * s%prad(k) / (s%peos(k) - s%prad(k))
            b_crit(k) = 4 * pi * s%rho(k) * s%csound(k) * s%csound(k) * Q * delta_grad
            b_crit(k) = b_crit(k) / (1 - Q * delta_grad + dgamma1(k))
            if (b_crit(k) < 0d0) b_crit(k) = 0d0
            b_crit(k) = sqrt(b_crit(k))
         end do

         names(1) = 'b_crit'
         do k = 1, nz
            vals(k,1) = b_crit(k)
         end do

      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr, profile_number=1
         real(dp) :: age
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         age = s% star_age

         if (age > profile_number*s% x_ctrl(1)) then
            s% need_to_save_profiles_now = .true.
            profile_number = profile_number + 1
         endif

         extras_finish_step = keep_going
      end function extras_finish_step


      subroutine my_other_energy(id, ierr)
         use const_def, only: pi
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i,j,k,l
         real(dp) :: true_irradiation, irradiation_dq, xq, log_irradiation, heat_efficiency
         real(dp) :: total_extra_counter, extra, qmin, qmax, qtop, qbot, q00, qp1

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Loop that finds deepest irradiated grid index location
         irradiation_dq = pi4*s% r(1)*s% r(1)*s% column_depth_for_irradiation/s% xmstar
         do i = 1, s% nz, 1
            xq = 1 - s% q(i)
            if (irradiation_dq < xq) exit
         end do

         ! Loop that finds 1 Mbar region location, where hydrogen metallic starts
         do j = 1, s% nz, 1
            if (s% peos(j) > 1e12) exit
         end do

         ! Loop that finds the layer of mass double that of index i (to get the luminosity in 
         ! a deeper region from the irradiated on, due to fluctuating values of luminosity)
         do k = 1, s% nz, 1
            if ((1 - s% q(k)) > 2*(1 - s% q(i))) exit
         end do
         
         ! Percentage that relates irradiation flux with extra deposition heat (from Thorngren & Fortney 2019)
         log_irradiation = log10( s%irradiation_flux / 1e9 )
         heat_efficiency = 2.37*exp( - ( log_irradiation - 0.14 )**2 / (.37)**2 / 2 )/100

         ! write(*, *) 'heat_efficiency', heat_efficiency

         ! Total deposited heat inside
         extra = heat_efficiency * s%irradiation_flux * pi * s%r(1) * s%r(1) / (s%m(j) - s%m(s%nz))


         ! Funtion that makes a progressive increase of the interior irradiation (to avoid relaxation)
         if (s%time_years < 1e6) extra = extra * sin(pi * s%time_years / 2 / 1e6 )

         ! write(*, *) 'time (years):', s%time_years
         
         total_extra_counter = 0

         qp1 = 0d0
         qmin = 0
         qmax = s%q(j)  ! This has to be the same j as the m(j) defined in extra
         do l= s%nz, 1, -1
            q00 = s%q(l)
            if (qp1 >= qmin .and. q00 <= qmax) then ! all inside of region
               s%extra_heat(l) = s% extra_heat(l) + extra
               total_extra_counter = total_extra_counter + extra * s%dm(l)
            else
               qtop = min(q00, qmax)
               qbot = max(qp1, qmin)
               if (qtop > qbot) then ! overlaps region
               write(*,*) "Hello"
                  s%extra_heat(l) = s%extra_heat(l) + extra * (qtop - qbot) / s%dq(l)
                  total_extra_counter = total_extra_counter + extra * s%dm(l) * (qtop - qbot) / s%dq(l)
               end if
            end if
            qp1 = q00
         end do

         !write(*, *) 'Energy deposited (erg s^-1): ', extra*(s%m(j) - s%m(s%nz)), total_extra_counter
         !write(*, '(a, F0.8)') 'Relative error (%): ', 100*abs(extra*(s%m(k) - s%m(s%nz)) - total_extra_counter)/extra/s% m(k)

      end subroutine my_other_energy

      end module run_star_extras
      
