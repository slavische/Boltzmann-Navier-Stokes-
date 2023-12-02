module cond_kin
   use constants
   use fluid_prop
   implicit none

contains

   function crit_radius(density, temperature) result(r_cr)
      real, intent(in) :: density, temperature
      real :: r_cr
      real :: s
      s = ss_ratio(density, temperature)
      if (s > 1.) then
         r_cr = 2.*surface_tension(temperature)*molar_mass/ &
                (rho_l*MOLAR_GAS_CONST*temperature*log(s))
      else
         r_cr = 0.
      end if
   end function crit_radius

   function nucleation_rate(density, temperature) result(j_cnt)
      implicit none
      real, intent(in) :: density, temperature
      real :: j_cnt

      double precision :: sigma, s, DeltaG, j0

      s = ss_ratio(density, temperature)
      if (s > 1.) then
         sigma = surface_tension(temperature)
         DeltaG = 16.*PI/3.*molar_volume**2.*sigma**3./ &
                  (BOLTZMANN_CONST*temperature*dlog(s))**2.
         DeltaG = 4./3.*sigma*PI*crit_radius(density, temperature)**2
         J0 = density**2/rho_l*dsqrt(2.*sigma/PI/molecule_mass/molecule_mass/molecule_mass)
         J_CNT = J0*dexp(-DeltaG/(BOLTZMANN_CONST*temperature))
      else
         J_CNT = 0.
      end if
   end function nucleation_rate

   function growth_rate(density, temperature) result(drdt)
      implicit none
      real, intent(in) :: density, temperature
      real :: drdt

      drdt = (pressure_IG_EOS(density, temperature) - sat_pressure(temperature))/rho_l/ &
             sqrt(2.*PI*MOLAR_GAS_CONST/molar_mass*temperature)
            if ( drdt < 0.) drdt = 0.
      !drdt = 0.
   end function growth_rate
end module cond_kin
