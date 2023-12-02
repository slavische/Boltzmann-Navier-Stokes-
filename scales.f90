module scales
   use fluid_prop
   implicit none
   save
!--SCALES------------------------------------
   real, private :: temperature_scale = 303., &
                    density_scale, &
                    velocity_scale, &
                    length_scale, &
                    time_scale
   logical, private :: scales_initialized = .false.
contains
   subroutine init_scales()
      if (.not. scales_initialized) then
         density_scale = sat_pressure(temperature_scale)/(gas_const*temperature_scale)
         velocity_scale = sqrt(gas_const*temperature_scale)
         length_scale = mean_free_path(density_scale)
         time_scale = length_scale/velocity_scale
         scales_initialized = .true.
      end if
   end subroutine init_scales
!--CONVERT TO PHYS---------------------------
   function phys_density(density_) result(density)
      real, intent(in) :: density_
      real :: density
      if (.not. scales_initialized) then
          call init_scales()
      end if
      density = density_*density_scale
   end function phys_density

   function phys_velocity(velocity_) result(velocity)
      real, intent(in) :: velocity_
      real :: velocity
      if (.not. scales_initialized) then
          call init_scales()
      end if
      velocity = velocity_*velocity_scale
   end function phys_velocity

   function phys_temperature(temperature_) result(temperature)
      real, intent(in) :: temperature_
      real :: temperature
      if (.not. scales_initialized) then
          call init_scales()
      end if
      temperature = temperature_*temperature_scale
   end function phys_temperature

   function phys_length(length_) result(length)
      real, intent(in) :: length_
      real :: length
      if (.not. scales_initialized) then
          call init_scales()
      end if
      length = length_*length_scale
   end function phys_length

   function phys_time(time_) result(time)
      real, intent(in) :: time_
      real :: time
      if (.not. scales_initialized) then
          call init_scales()
      end if
      time = time_*time_scale
   end function phys_time

!--COVERT TO DIMENSIONLESS-------------------
   function dl_density(density) result(density_)
      real, intent(in) :: density
      real :: density_
      if (.not. scales_initialized) then
          call init_scales()
      end if
      density_ = density/density_scale
   end function dl_density

   function dl_velocity(velocity) result(velocity_)
      real, intent(in) :: velocity
      real :: velocity_
      if (.not. scales_initialized) then
          call init_scales()
      end if
      velocity_ = velocity/velocity_scale
   end function dl_velocity

   function dl_temperature(temperature) result(temperature_)
      real, intent(in) :: temperature
      real :: temperature_
      if (.not. scales_initialized) then
          call init_scales()
      end if
      temperature_ = temperature/temperature_scale
   end function dl_temperature

   function dl_length(length) result(length_)
      real, intent(in) :: length
      real :: length_
      
      if (.not. scales_initialized) then
          call init_scales()
      end if
      
      length_ = length/length_scale
   end function dl_length

   function dl_time(time) result(time_)
      real, intent(in) :: time
      real :: time_
      if (.not. scales_initialized) then
          call init_scales()
      end if
      time_ = time/time_scale
   end function dl_time

end module scales
