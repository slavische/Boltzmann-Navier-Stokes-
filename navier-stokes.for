	module navier_stokes
      use utils
      implicit none
      contains
      
      subroutine flowns(rl_b,tl_b,ul_b,rhoa,temp,vel,cva,n)
	integer, parameter :: rx=5002
	integer i,n,ktt
	real gama,cva,umix0,tmix0,rhomix0
	real rhoa(rx),temp(rx),vel(rx)
	real rl_b,tl_b,ul_b
	
      n=960

	ktt=1
	umix0=0.
	tmix0=1.
	rhomix0=1.


c ***** <a-Ar> **** <b-He> ******
c	mamb=40./4.
c	dadb=3.67/2.2
c
	gama=9./7.
c	gamb=5./3.
c ***** <a-Ar> **** <b-He> ******
c ***** <a-H20> ***** <b-N2> *****
c	mamb=2.99/4.65
c	dadb=4.34/3.7

c	gama=9./7.
c	gamb=7./5.
c ***** <a-H20> ***** <b-N2> *****
	cva=1./(gama-1.)


          IF(Ktt.EQ.2) THEN
C******************************************************
             OPEN(15,FILE='init1_',STATUS='OLD',
     *       ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
             READ(15) rhoa,temp,vel
             CLOSE(15)
C******************************************************
		else
	do i=1, n
		rhoa(i)=rl_b
		vel(i)=ul_b
		temp(i)=tl_b
	enddo
          END IF


      
	open(1,file='init.dat',status='unknown')

	do i=1,n,1
 
      write(1,99) rhoa(i),temp(i),vel(i)
   99            FORMAT(2F9.5,f11.7)

	end do

	close(2)

	return
	end

	subroutine fprint_ns(file_num,xx,rhoa,temp,vel,n)
		integer, parameter :: rx=5002
c	parameter (tend=3.47,kk=100)
		character(len=10), intent(in) :: file_num
		real, intent(in) :: xx(rx), rhoa(rx), temp(rx), vel(rx)
		integer, intent(in) :: n
		integer :: j

		character(len=30) :: file_name
		file_name = 'ev-con-'//trim(adjustl(file_num))//'.dat'
		open(2, file=trim(file_name))

		write(2, '(5(A12, 2x))') 'xx/10.', 'rhoa', 'temp', 'vel', 'vel*rhoa'
		write(2, '(5(A12, 2x))') '[-]', '[-]','[-]', '[-]','[-]'
		do j=1,n,1
			write(2,'(5(Es12.5, 2x))')xx(j)/10.,
     *	rhoa(j),temp(j),vel(j),vel(j)*rhoa(j)

		end do

	close(2)

	return
	end


	subroutine navie_stoks(rl_b,tl_b,ul_b,rl_ns,tl_ns,ul_ns,
     *	dtau,rhoa,temp,cva,vel,xx,n,rr_b,tr_b,ur_b,rr_ns,tr_ns,ur_ns)


	real cva,cpa,ca,dx,time,dtau,umix0,tmix
	real rl_b,tl_b,ul_b,rl_ns,tl_ns,ul_ns
	real rr_b,tr_b,ur_b,rr_ns,tr_ns,ur_ns
c	parameter (rx=1002,pi=3.1415926)
	integer, parameter :: rx=5002
c	parameter (tend=3.47,kk=100)
	integer j,n,i,jr,usl,usp
	real kk1,kk2
	real ei,gama,vsound,vmax
	real sumrha,pi,sumrhoe
	real rhoa(rx),rhou(rx),rhoe(rx),temp(rx),vel(rx),xx(rx)
	real rho(rx),preje(rx),rho1(rx)
	real a(rx),b(rx),l(rx)
	real a1(rx),c1(rx),b1(rx),f1(rx),y(0:rx),dxx(rx)
	real l1(rx),km(rx),kpl(rx),d1(rx),diffm(rx)
	real ulef,urig,wlef,wrig,blef,brig
	real rhol,rhor,tlef,trig,cmix(rx),cvmd(rx),cvmix(rx),kmd(rx)
	real rhoal,rhoar,rhoalef,rhoarig
	real tl,tp,ulefa,uriga,caa(rx)
	real flowdiffa(rx),kmdif(rx),kmu(rx),kmp(rx)
	real dudx(rx),flowdiffa1(rx)
	real kmda(rx),d2rhadx(rx)

c	kk1=2.
c	kk2=700.
	n=960

c	rhoalef=1.974
c	rhoarig=1.026

c	rhoalef=1.1627
c	rhoarig=0.995

	tlef=tl_b
	trig=tr_b


c      open(2,file='he_tepl.dat',status='unknown')
c	open(4,file='voda-azotflow.dat',status='unknown')  

	pi=3.1415926

	jr=1

c	umix0=0.

c	dx=0.165
	dx=0.5
c	dx=0.3

c	tmix=1.

c	dtau=0.03
c	dtau=0.0001

c	mamb=10./1.
c	dadb=3.67/2.2

c ***** <a-Ar> **** <b-He> ******
c	mamb=40./4.
c	dadb=3.67/2.2
c
c	gama=5./3.
c	gamb=5./3.
c ***** <a-Ar> **** <b-He> ******
c ***** <a-H20> ***** <b-N2> *****
c	mamb=2.99/4.65
c	dadb=4.34/3.7

c	gama=9./7.
c	gamb=7./5.
c ***** <a-H20> ***** <b-N2> *****
c	cva=1./(gama-1.)
c	cvb=mamb/(gamb-1.)


c	vsound=dsqrt(gamb*tmix)
c	dtau=0.3*dx/dabs(vsound)/kk1

	ulefa=ul_b
	uriga=ur_b

	xx(1)=10.+dx/2.
	do j=2,n
	xx(j)=xx(j-1)+dx
	enddo

c	dxx(1)=dx/2.

	do j=1,n+1
	dxx(j)=dx
	enddo

c	dxx(n+1)=dx/2.

	
	do i=1, n
c		rhoa(i)=1.0
c		u(i)=umix0
		rhou(i)=rhoa(i)*vel(i)
c		temp(i)=tmix
		rhoe(i)=cva*temp(i)*rhoa(i)
	preje(i)=rhoa(i)*temp(i)
	enddo


c	time=0.



	rhoal=rl_b
	rhoar=rr_b

c
		do j=1,n
		a(j)=rhoe(j)
		b(j)=vel(j)
		l(j)=rhoa(j)*temp(j)
		enddo

		blef=ul_b
		brig=ur_b

c		wlef=rhoal*tlef*cva*ulefa*dtau/dx
		wlef=rl_b*tl_b*cva*ul_b*dtau/dx
		wrig=rr_b*tr_b*cva*ur_b*dtau/dx
c		wrig=rhoar*trig*cva*uriga*dtau/dx
		
		call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)
	
	do j=1,n
	rhoe(j)=a(j)
	enddo

c ****************************************************
		do j=1,n
		a(j)=rhou(j)
		b(j)=rhoa(j)*temp(j)
		l(j)=1.
		enddo


		blef=rl_b*tl_b
		brig=rr_b*tr_b


c		wlef=rhoal*ulefa*ulefa*dtau/dx
		wlef=rl_b*ul_b*ul_b*dtau/dx
		wrig=rr_b*ur_b*ur_b*dtau/dx
c		wrig=rhoar*uriga*uriga*dtau/dx


			call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)
	
	do j=1,n
	rhou(j)=a(j)
	enddo

c ******************************************************	

		do j=1,n
		a(j)=rhoa(j)
		b(j)=0.
		l(j)=1.
		enddo
		
		blef=0.	
		brig=0.
		

c		wlef=rhoal*ulefa*dtau/dx
		wlef=rl_b*ul_b*dtau/dx
		wrig=rr_b*ur_b*dtau/dx
c		wrig=rhoar*uriga*dtau/dx
	
			call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)

		do j=1,n
		rhoa(j)=a(j)
		enddo

c ******************************************************	
		do j=1,n
		vel(j)=rhou(j)/rhoa(j)
		temp(j)=rhoe(j)/rhoa(j)/cva
		enddo

c *****************momentum equation **************************
cc
	do j=2,n
	kmu(j)=5./16.*sqrt(2.*pi)*(sqrt(temp(j))+
     *	sqrt(temp(j-1)))/2.
	enddo

	do j=1,n
	l1(j)=rhoa(j)
	y(j)=vel(j)
	d1(j)=0.
	enddo
	
	usl=1
	usp=1

c	tl=2.*u(1)-u(2)
c	tp=2.*u(n)-u(n-1)

c	tl=ulef
c	tp=urig

	tl=ul_b
	tp=ur_b
	
	dtau=4./3.*dtau
c
		call teplo(kmu,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)

	do j=1,n
	vel(j)=y(j)
	enddo

c	rl_ns=rhoa(1)
c	ul_ns=y(0)

	dtau=3./4.*dtau
c	kmu=km

c **********energy equation ************************************	
	usl=1
	usp=1

	tl=tl_b
	tp=tr_b
	
	do j=2,n
	kmp(j)=5./16.*sqrt(2.*pi)*5./2.*cva*(sqrt(temp(j))+
     *	sqrt(temp(j-1)))/2.

	enddo

	do j=1,n
	l1(j)=rhoa(j)*cva
	y(j)=temp(j)
	d1(j)=0.
	enddo

	do j=2,n-1
	dudx(j)=(vel(j+1)-vel(j-1))/2./dxx(j)
	enddo

	dudx(1)=(vel(2)-ul_b)/2./dxx(2)
	dudx(n)=(ur_b-vel(n-1))/2./dxx(n-1)

	call teplo(kmp,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)

	do j=1,n
	temp(j)=y(j)
	enddo
	
	do j=1,n
	temp(j)=dtau*(4./3.*kmu(j)*dudx(j)*dudx(j))+temp(j)
	enddo

c **************************************************


c	do i=1, n
c		rhou(i)=rhoa(i)*u(i)
c		rhoe(i)=rhoa(i)*temp(i)*cva
c	enddo

	rl_ns=rhoa(1)
	tl_ns=temp(1)
	ul_ns=vel(1)

	rr_ns=rhoa(n)
	tr_ns=temp(n)
	ur_ns=vel(n)

	return
      end

	subroutine teplo(km,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)
	integer, parameter :: rx=5002
	real a1(rx),b1(rx),c1(rx),l1(rx),alpha(rx),betta(rx)
	real y(0:rx),f1(rx),d1(rx),km(rx),dxx(rx)
	real dtau,tl,tp
	real kappa2,kappa1
	integer j,n,usl,usp

	
	km(1)=km(2)
c	km(n)=km(n-1)
	km(n+1)=km(n)
	do j=1,n
	a1(j)=2.*dtau*km(j)/dxx(j)/(dxx(j+1)+dxx(j))
	b1(j)=2.*dtau*km(j+1)/dxx(j+1)/(dxx(j+1)+dxx(j))
	c1(j)=a1(j)+b1(j)+l1(j)
	f1(j)=l1(j)*y(j)+d1(j)
	enddo

	if(usl==2)then
	alpha(1)=1.
	betta(1)=0.
	else
	alpha(1)=0.
	betta(1)=tl
	endif

	do j=1,n
	alpha(j+1)=b1(j)/(c1(j)-a1(j)*alpha(j))
	betta(j+1)=(a1(j)*betta(j)+f1(j))/(c1(j)-a1(j)*alpha(j))
	enddo

	if(usp==2)then
	kappa2=1.
	kappa1=0.
	else
	kappa2=0.
	kappa1=tp
	endif

	y(n+1)=(kappa2*betta(n+1)+kappa1)/(1.-kappa2*alpha(n+1))
	do j=n,0,-1
	y(j)=alpha(j+1)*y(j+1)+betta(j+1)
	enddo



	return
      end subroutine teplo
      
      
      subroutine navie_stoks_cond(n, dtau,
     *                            rl_b, tl_b, ul_b, rr_b, tr_b, ur_b,
     *                            rhoa, temp, vel, cva,
     *                            Omega0, Omega1, Omega2, C_l,
     *                            xx, rl_ns,tl_ns,ul_ns, 
     *                            rr_ns,tr_ns,ur_ns)
      use scales
	  use mom_eq
      implicit none
      integer, intent(in) :: n
      real, intent(in) :: rl_b,tl_b,ul_b
      real, intent(in) :: rr_b,tr_b,ur_b
      real, intent(in) :: cva
      real, intent(inout) :: dtau
      real, dimension(:), intent(inout) :: rhoa, temp, vel
      real, dimension(:), intent(inout) :: Omega0, Omega1, Omega2, C_l
      real, dimension(:), intent(out) ::  xx
      real, intent(out) :: rl_ns, tl_ns, ul_ns
      real, intent(out) :: rr_ns, tr_ns, ur_ns
      
	real dx
      real, dimension(:), allocatable ::  rho, rhou, rhoe, preje, 
     *                                    rhoOmega0, rhoOmega1, 
     *                                    rhoOmega2, rhoC_l
	real, dimension(:), allocatable :: a, b, l
      real, dimension(:), allocatable :: rhoaPhys, rhoPhys, velPhys, 
     *                                   tempPhys
      real :: dxPhys, dtauPhys, ulPhys, urPhys
      integer j,i,jr,usl,usp
      integer, parameter :: rx=5002
	real a1(rx),c1(rx),b1(rx),f1(rx),y(0:rx),dxx(rx)
	real l1(rx),km(rx),kpl(rx),d1(rx),diffm(rx)
	real ulef,urig,wlef,wrig,blef,brig
	real rhor,tlef,trig,cmix(rx),cvmd(rx),cvmix(rx),kmd(rx)
	real rhoal,rhoar,rhoalef,rhoarig
	real tl,tp,ulefa,uriga,caa(rx)
	real flowdiffa(rx),kmdif(rx),kmu(rx),kmp(rx)
	real dudx(rx),flowdiffa1(rx)
	real kmda(rx),d2rhadx(rx)

      allocate(rho(n), rhou(n), rhoe(n), preje(n), 
     *         rhoOmega0(n), rhoOmega1(n), rhoOmega2(n), rhoC_l(n))
      allocate(a(n), b(n), l(n))
      allocate(rhoaPhys(n), rhoPhys(n), velPhys(n), tempPhys(n))
      
	tlef=tl_b
	trig=tr_b

	jr=1

	ulefa=ul_b
	uriga=ur_b
      
      dx=0.5
	xx(1)=10.+dx/2.
	do j=2,n
	    xx(j)=xx(j-1)+dx
	end do

	do j=1,n+1
	    dxx(j)=dx
	end do

	do i=1, n
          rho(i) = rhoa(i)/(1. - C_l(i))
		rhou(i) = rho(i)*vel(i)
		rhoe(i) = cva*temp(i)*rho(i)
	    preje(i) = rho(i)*temp(i)
	end do

	rhoal=rl_b
	rhoar=rr_b

      do j=1,n
		a(j)=rhoe(j)
		b(j)=vel(j)
		l(j)=rho(j)*temp(j)
      enddo

      blef=ul_b
      brig=ur_b

      wlef=rl_b*tl_b*cva*ul_b*dtau/dx
      wrig=rho(n)*tr_b*cva*ur_b*dtau/dx
		
      call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)
	
	do j=1,n
	    rhoe(j)=a(j)
	end do

c ****************************************************
		do j=1,n
		a(j)=rhou(j)
		b(j)=rho(j)*temp(j)
		l(j)=1.
		enddo


		blef=rl_b*tl_b
		brig=rr_b*tr_b


		wlef=rl_b*ul_b*ul_b*dtau/dx
		wrig=rho(n)*ur_b*ur_b*dtau/dx

          call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)
	
	do j=1,n
	rhou(j)=a(j)
      enddo
      
c     MOMENT EQ
      dtauPhys = phys_time(dtau)
      dxPhys = phys_length(dx)
      ulPhys = phys_velocity(ul_b)
      urPhys = phys_velocity(ur_b)
      do i=1, n
          rhoPhys(i) = phys_density(rhoa(i))/(1. - C_l(i))
          velPhys(i) = phys_velocity(vel(i))
          rhoOmega0(i) = rhoPhys(i)*Omega0(i)
          rhoOmega1(i) = rhoPhys(i)*Omega1(i)
          rhoOmega2(i) = rhoPhys(i)*Omega2(i)
          rhoC_l(i) = rhoPhys(i)*C_l(i)
      end do
      
      call mom_transport(n, dxPhys, dtauPhys, ulPhys, urPhys, velPhys,
     *rhoOmega0, rhoOmega1, rhoOmega2, rhoC_l)

      

c ******************************************************	

		do j=1,n
		a(j)=rho(j)
		b(j)=0.
		l(j)=1.
		enddo
		
		blef=0.	
		brig=0.
		

		wlef=rl_b*ul_b*dtau/dx
		wrig=rho(n)*ur_b*dtau/dx
	
          call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)

      do j=1,n
		rho(j)=a(j)
      enddo

c ******************************************************	
		do j=1,n
              vel(j)=rhou(j)/rho(j)
              temp(j)=rhoe(j)/rho(j)/cva
              Omega0(j) = rhoOmega0(j)/phys_density(rho(j))
              Omega1(j) = rhoOmega1(j)/phys_density(rho(j))
              Omega2(j) = rhoOmega2(j)/phys_density(rho(j))
              C_l(j) = rhoC_l(j)/phys_density(rho(j))
			  rhoa(j) = rho(j)*(1. -C_l(j))
		enddo

c *****************momentum equation **************************
cc
	do j=2,n
	kmu(j)=5./16.*sqrt(2.*pi)*(sqrt(temp(j))+
     *	sqrt(temp(j-1)))/2.
	enddo

	do j=1,n
	l1(j)=rho(j)
	y(j)=vel(j)
	d1(j)=0.
	enddo
	
	usl=1
	usp=1

c	tl=2.*u(1)-u(2)
c	tp=2.*u(n)-u(n-1)

c	tl=ulef
c	tp=urig

	tl=ul_b
	tp=ur_b
	
	dtau=4./3.*dtau
c
        
		call teplo(kmu,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)

	do j=1,n
	vel(j)=y(j)
	enddo

c	rl_ns=rhoa(1)
c	ul_ns=y(0)

	dtau=3./4.*dtau
c	kmu=km

c **********energy equation ************************************	
	usl=1
	usp=1

	tl=tl_b
	tp=tr_b
	
	do j=2,n
	kmp(j)=5./16.*sqrt(2.*pi)*5./2.*cva*(sqrt(temp(j))+
     *	sqrt(temp(j-1)))/2.

	enddo

	do j=1,n
	l1(j)=rho(j)*cva
	y(j)=temp(j)
	d1(j)=0.
	enddo

	do j=2,n-1
	dudx(j)=(vel(j+1)-vel(j-1))/2./dxx(j)
	enddo

	dudx(1)=(vel(2)-ul_b)/2./dxx(2)
	dudx(n)=(ur_b-vel(n-1))/2./dxx(n-1)

	call teplo(kmp,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)

	do j=1,n
	temp(j)=y(j)
	enddo
	
	do j=1,n
	temp(j)=dtau*(4./3.*kmu(j)*dudx(j)*dudx(j))+temp(j)
	enddo

c **************************************************


c	do i=1, n
c		rhou(i)=rhoa(i)*u(i)
c		rhoe(i)=rhoa(i)*temp(i)*cva
c	enddo


	  do i=1, n
      rhoaPhys(i) = phys_density(rhoa(i))
      tempPhys(i) = phys_temperature(temp(i))
	  end do
      call sources(n, dtauPhys, rhoaPhys, tempPhys, 
     *omega0, omega1, omega2, C_l)
	  do i=1, n
	  rhoa(i) = dl_density(rhoaPhys(i))
	  temp(i) = dl_temperature(tempPhys(i))
	  end do



	rl_ns=rhoa(1)
	tl_ns=temp(1)
	ul_ns=vel(1)

	rr_ns=rhoa(n)
	tr_ns=temp(n)
	ur_ns=vel(n)

	return
      end subroutine navie_stoks_cond

	  subroutine sources(n, dtau, rhoa, temp, omega0, omega1, 
     *omega2, C_l)
		use constants
		use scales

		implicit none
		integer, intent(in) :: n
		! PHYS
	  real, intent(in) :: dtau
      real, dimension(:), intent(inout) :: rhoa, temp, omega0, 
     *omega1, omega2, C_l

		real, dimension(:), allocatable :: rho
		real, dimension(:), allocatable :: s0, s1, s2, s3
		integer :: j
		real :: J_CNT, drdt, r_cr, dtau_, tmp, tau
  
		real,parameter :: EPS = 1.e-06


		allocate (rho(n))
		allocate (s0(n), s1(n), s2(n), s3(n))


			call rhs(n, rhoa, temp, omega0, omega1, omega2, 
     *s0, s1, s2, s3)

			do j = 1, n
				rho(j) = rhoa(j)/(1.-C_l(j))
				omega0(j) = omega0(j) + dtau/rho(j)*s0(j)
				omega1(j) = omega1(j) + dtau/rho(j)*s1(j)
				omega2(j) = omega2(j) + dtau/rho(j)*s2(j)
				c_l(j) = c_l(j) + dtau/rho(j)*s3(j)
				temp(j) = temp(j) + dtau*heat_evp/rho(j)/Cv_v*s3(j)
				rhoa(j) =rho(j)*(1.-C_l(j))
			end do

	 end subroutine sources

	 	subroutine rhs(n, density, temperature, omega0, omega1, omega2,
     *s0, s1, s2, s3)
	 		use constants
			use cond_kin
			implicit none
			integer, intent(in) :: n
			real, dimension(:), intent(in) :: density, temperature
			real, dimension(:), intent(in) :: omega0, omega1, omega2
			real, dimension(:), intent(out) :: s0, s1, s2, s3
			integer :: j
			real :: J_CNT, drdt, r_cr

			do j = 1, n
				J_CNT = nucleation_rate(density(j), temperature(j))
				drdt = growth_rate(density(j), temperature(j))
				r_cr = crit_radius(density(j), temperature(j))
				s0(j) = J_CNT 
				s1(j) = J_CNT*r_cr + drdt*omega0(j)
				s2(j) = J_CNT*r_cr**2. + 2.*drdt*omega1(j)
				s3(j) = J_CNT*r_cr**3. + 3.*drdt*omega2(j)
				s3(j) = 4./3.*PI*rho_l*s3(j)
		 	end do

		end subroutine rhs


      	subroutine print_ss(file_num, n, xx, rhoa, temp)
			use scales
			use fluid_prop
			implicit none
			character(len=10), intent(in) :: file_num
			integer, intent(in) :: n
			real, dimension(:), intent(in) :: xx, rhoa, temp
			integer :: i 

			character(len=30) :: file_name
			file_name = 'ss-ratio-'//trim(adjustl(file_num))//'.dat'
			open(1, file=trim(file_name))

            write(1, '(2(A12, 2x))') 'xx/10.', 'ss_ratio'
            write(1, '(2(A12, 2x))') '[-]', '[-]'
            do i = 1, n
                write(1, '(2(Es12.5, 2x))') xx(i)/10., 
     *ss_ratio(phys_density(rhoa(i)), phys_temperature(temp(i)))
            end do
			
            close(1)
		end subroutine print_ss
	end module navier_stokes