module fluid_prop
   use constants
   implicit none

   real, parameter :: molar_mass = 18.0153e-03, &
                      cross_section = PI*3.4e-10**2, &
                      rho_l = 995.65, &
                      Cp_l = 4180.1, &
                      Cp_v = 1917.8, &
                      Cv_v = 1445.1, &
                      heat_evp = 2430.2e+03
!  SATURATION LINE
   real, parameter :: AA = 4.6543, BB = 1435.264, CC = -64.848
   real, parameter :: gas_const = MOLAR_GAS_CONST/molar_mass, &
                      molecule_mass = molar_mass/AVOGADRO_CONST, &
                      molar_volume = molar_mass/rho_l
!  SIGMA DATA
   integer, parameter :: size_data = 73
   real, dimension(size_data), parameter :: temp_data = [280.00, &
                                                         285.00, &
                                                         290.00, &
                                                         295.00, &
                                                         300.00, &
                                                         305.00, &
                                                         310.00, &
                                                         315.00, &
                                                         320.00, &
                                                         325.00, &
                                                         330.00, &
                                                         335.00, &
                                                         340.00, &
                                                         345.00, &
                                                         350.00, &
                                                         355.00, &
                                                         360.00, &
                                                         365.00, &
                                                         370.00, &
                                                         375.00, &
                                                         380.00, &
                                                         385.00, &
                                                         390.00, &
                                                         395.00, &
                                                         400.00, &
                                                         405.00, &
                                                         410.00, &
                                                         415.00, &
                                                         420.00, &
                                                         425.00, &
                                                         430.00, &
                                                         435.00, &
                                                         440.00, &
                                                         445.00, &
                                                         450.00, &
                                                         455.00, &
                                                         460.00, &
                                                         465.00, &
                                                         470.00, &
                                                         475.00, &
                                                         480.00, &
                                                         485.00, &
                                                         490.00, &
                                                         495.00, &
                                                         500.00, &
                                                         505.00, &
                                                         510.00, &
                                                         515.00, &
                                                         520.00, &
                                                         525.00, &
                                                         530.00, &
                                                         535.00, &
                                                         540.00, &
                                                         545.00, &
                                                         550.00, &
                                                         555.00, &
                                                         560.00, &
                                                         565.00, &
                                                         570.00, &
                                                         575.00, &
                                                         580.00, &
                                                         585.00, &
                                                         590.00, &
                                                         595.00, &
                                                         600.00, &
                                                         605.00, &
                                                         610.00, &
                                                         615.00, &
                                                         620.00, &
                                                         625.00, &
                                                         630.00, &
                                                         635.00, &
                                                         640.00]
   real, dimension(size_data), parameter :: surface_tension_data = [0.074677, &
                                                                    0.073951, &
                                                                    0.073210, &
                                                                    0.072455, &
                                                                    0.071686, &
                                                                    0.070903, &
                                                                    0.070106, &
                                                                    0.069295, &
                                                                    0.068470, &
                                                                    0.067632, &
                                                                    0.066781, &
                                                                    0.065917, &
                                                                    0.065040, &
                                                                    0.064150, &
                                                                    0.063248, &
                                                                    0.062333, &
                                                                    0.061406, &
                                                                    0.060467, &
                                                                    0.059517, &
                                                                    0.058555, &
                                                                    0.057581, &
                                                                    0.056596, &
                                                                    0.055601, &
                                                                    0.054595, &
                                                                    0.053578, &
                                                                    0.052551, &
                                                                    0.051514, &
                                                                    0.050468, &
                                                                    0.049411, &
                                                                    0.048346, &
                                                                    0.047272, &
                                                                    0.046189, &
                                                                    0.045098, &
                                                                    0.043999, &
                                                                    0.042891, &
                                                                    0.041777, &
                                                                    0.040655, &
                                                                    0.039527, &
                                                                    0.038392, &
                                                                    0.037252, &
                                                                    0.036105, &
                                                                    0.034954, &
                                                                    0.033797, &
                                                                    0.032637, &
                                                                    0.031472, &
                                                                    0.030304, &
                                                                    0.029133, &
                                                                    0.027959, &
                                                                    0.026784, &
                                                                    0.025608, &
                                                                    0.024430, &
                                                                    0.023253, &
                                                                    0.022077, &
                                                                    0.020902, &
                                                                    0.019730, &
                                                                    0.018561, &
                                                                    0.017396, &
                                                                    0.016236, &
                                                                    0.015082, &
                                                                    0.013937, &
                                                                    0.012800, &
                                                                    0.011673, &
                                                                    0.010559, &
                                                                    0.0094591, &
                                                                    0.0083756, &
                                                                    0.0073112, &
                                                                    0.0062690, &
                                                                    0.0052528, &
                                                                    0.0042676, &
                                                                    0.0033194, &
                                                                    0.0024169, &
                                                                    0.0015728, &
                                                                    0.00080882]

contains

   function pressure_IG_EOS(density, temperature) result(p)
      implicit none
      real, intent(in) :: density, temperature
      real :: p

      p = density*MOLAR_GAS_CONST/molar_mass*temperature
   end function pressure_IG_EOS

   function sat_pressure(temperature) result(p_s)
      real, intent(in) :: temperature
      real :: p_s

      p_s = 1.e+5*10.**(AA - (BB/(temperature + CC)))
   end function sat_pressure

   function ss_ratio(density, temperature) result(s)
      real, intent(in) :: density, temperature
      real :: s

      s = pressure_IG_EOS(density, temperature)/sat_pressure(temperature)
   end function ss_ratio

   function surface_tension(temperature) result(sigma)
      real, intent(in) :: temperature
      real :: sigma

      sigma = linear_interpolation(temp_data, surface_tension_data, temperature)
   end function surface_tension

   real function mean_free_path(density) result(path_length)
      real, intent(in) :: density

      path_length = 1.0/(cross_section*density/molar_mass*AVOGADRO_CONST)
   end function mean_free_path

   function linear_interpolation(xData, yData, xVal) result(yVal)
      real, dimension(:), intent(in) :: xData, yData
      real, intent(in) :: xVal
      real :: yVal
      integer :: i, n

      n = size(xData)

      if (xVal <= xData(1)) then
         yVal = yData(1)
      else if (xVal >= xData(n)) then
         yVal = yData(n)
      else
         do i = 2, n
            if (xVal <= xData(i)) then
               yVal = yData(i - 1) + (yData(i) - yData(i - 1))*(xVal - xData(i - 1))/(xData(i) - xData(i - 1))
               exit
            end if
         end do
      end if
   end function linear_interpolation

end module fluid_prop
