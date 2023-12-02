
    module mom_eq
        implicit none

    contains
    
        subroutine flow_mom(n, Omega0, Omega1, Omega2, C_l)
        use constants
        use fluid_prop
			implicit none
			integer, intent(in) :: n
            real, dimension(:), intent(inout) :: Omega0, Omega1, Omega2, C_l
			integer :: i
        
            do i = 1, n
                Omega0(i) = 0.
                Omega1(i) = 0.
                Omega2(i) =0.
                C_l(i) = 0.
            end do
            
            open(1,file='mom_init.dat',status='unknown')
            write(1, '(5(A12, 2x))') 'Omega0','Omega1', 'Omega2', 'C_l'
            write(1, '(5(A12, 2x))') '[kg^-1]','[m kg^-1]', '[m^2 kg^-1]', '[-]'
            do i = 1, n
                write(1, '(4(Es12.5, 2x))') Omega0(i), Omega1(i), Omega2(i), C_l(i)
            end do
            close(1)
        end subroutine flow_mom
       
        subroutine fprint_mom(file_num, n, xx, Omega0, Omega1, Omega2, C_l)
            implicit none
            character(len=10), intent(in) :: file_num
            real, dimension(:), intent(in) :: xx
            integer, intent(in) :: n
            real, dimension(:), intent(in) :: Omega0, Omega1, Omega2, C_l

            integer :: i
			character(len=30) :: file_name
			file_name = 'mom-'//trim(adjustl(file_num))//'.dat'
			open(1, file=trim(file_name))

            write(1, '(6(A12, 2x))') 'xx/10.', 'Omega0','Omega1', 'Omega2', 'C_l', 'r_avg'
            write(1, '(6(A12, 2x))') '[-]', '[kg^-1]','[m kg^-1]', '[m^2 kg^-1]', '[-]', '[m]'
            do i = 1, n
                write(1, '(6(Es12.5, 2x))') xx(i)/10., Omega0(i), Omega1(i), Omega2(i), C_l(i), r_avg(Omega0(i), C_l(i))
            end do
            close(1)
        
        end subroutine fprint_mom
        
        subroutine mom_transport(n, dx, dtau, ul, ur, vel, &
                                    rhoOmega0, rhoOmega1, rhoOmega2, rhoC_l)
        use utils    
        implicit none
            integer, intent(in) :: n
            ! PHYS
            real, intent(in) :: dtau, dx
            real, intent(in) :: ul
	        real, intent(in) :: ur
            real, dimension(:), intent(in) :: vel
            real, dimension(:), intent(inout) :: rhoOmega0, rhoOmega1, rhoOmega2, rhoC_l
            real, dimension(:), allocatable :: a, b, l
            real :: blef, brig, wlef, wrig
            integer :: j

            allocate (a(n), b(n), l(n))

        !   LIQUID MASS FRACTION EQUATION
            do j = 1, n
                a(j) = rhoC_l(j)
                b(j) = 0.
                l(j) = 1.
            end do

            blef = 0.
            brig = 0.

            wlef = 0.
            wrig = rhoC_l(n)*ur*dtau/dx

            call shasta(a, b, l, dx, vel, dtau, n, blef, brig, wlef, wrig)

            do j = 1, n
                rhoC_l(j) = a(j)
            end do

        !   2 MOMENT EQUATION
            do j = 1, n
                a(j) = rhoOmega2(j)
                b(j) = 0.
                l(j) = 1.
            end do

            blef = 0.
            brig = 0.

            wlef = 0.
            wrig = rhoOmega2(n)*ur*dtau/dx

            call shasta(a, b, l, dx, vel, dtau, n, blef, brig, wlef, wrig)

            do j = 1, n
                rhoOmega2(j) = a(j)
            end do

        !   1 MOMENT EQUATION
            do j = 1, n
                a(j) = rhoOmega1(j)
                b(j) = 0.
                l(j) = 1.
            end do

            blef = 0.
            brig = 0.

            wlef = 0.
            wrig = rhoOmega1(n)*ur*dtau/dx

            call shasta(a, b, l, dx, vel, dtau, n, blef, brig, wlef, wrig)

            do j = 1, n
                rhoOmega1(j) = a(j)
            end do

        !   0 MOMENT EQUATION
            do j = 1, n
                a(j) = rhoOmega0(j)
                b(j) = 0.
                l(j) = 1.
            end do

            blef = 0.
            brig = 0.

            wlef = ul
            wrig = rhoOmega0(n)*ur*dtau/dx

            call shasta(a, b, l, dx, vel, dtau, n, blef, brig, wlef, wrig)

            do j = 1, n
                rhoOmega0(j) = a(j)
            end do
            
        end subroutine mom_transport
        
        function r_avg(Omega0, C_l) result(r)
            use constants
            use fluid_prop
            real, intent(in) :: Omega0, C_l
            real :: r
        
            r = (3.*C_l/(4.*PI*rho_l*Omega0))**(1./3.)
        end function r_avg
    end module mom_eq