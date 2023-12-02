      module utils
      implicit none
      contains
      
      subroutine shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)
	integer j,n
	integer, parameter :: rx=5002
	real dtau,dtdx,dx
	real qplus(rx),qmin(rx),splus(rx),smin(rx),daudx(rx),dbdx(rx)
	real a1(rx),a2(rx),a3(rx),a4(rx)
	real a(rx),b(rx),l(rx),vel(rx)
	real ac(rx),dp(rx),dpl(rx),fc(rx)
	real koef1,koef2,koef3,koef4,koef5
	real blef,brig,wlef,wrig
	
	dtdx=dtau/dx
	do j=2,n-1


	qplus(j)=(1./2.-vel(j)*dtdx)/(1.+(vel(j+1)-vel(j))*dtdx)
	qmin(j)=(1./2.+vel(j)*dtdx)/(1.-(vel(j-1)-vel(j))*dtdx)

	splus(j)=l(j)*(b(j+1)-b(j))*dtdx
	smin(j)=l(j)*(b(j)-b(j-1))*dtdx

	a1(j)=qplus(j)*qplus(j)*(a(j+1)-a(j))/2.
	a2(j)=qmin(j)*qmin(j)*(a(j-1)-a(j))/2.
	a3(j)=qplus(j)*(a(j)-splus(j))
	a4(j)=qmin(j)*(a(j)-smin(j))
	enddo

	do j=2,n-1
		ac(j)=a1(j)+a2(j)+a3(j)+a4(j)
c		a(j)=a1(j)+a2(j)+a3(j)+a4(j)
	enddo



	qplus(1)=(1./2.-vel(1)*dtdx)/(1.+(vel(2)-vel(1))*dtdx)
	splus(1)=l(1)*(b(2)-blef)*dtdx
	a1(1)=qplus(1)*qplus(1)*(a(2)-a(1))/2.
	a3(1)=qplus(1)*(a(1)-splus(1))
	
c	a(1)=a1(1)+a3(1)+a(1)/2.+wlef
	a(1)=a1(1)+a3(1)+a(1)/2.+wlef+
     * (blef-b(1))*l(1)*dtdx

	qmin(n)=(1./2.+vel(n)*dtdx)/(1.-(vel(n-1)-vel(n))*dtdx)
	smin(n)=l(n)*(brig-b(n-1))*dtdx
	a2(n)=qmin(n)*qmin(n)*(a(n-1)-a(n))/2.
	a4(n)=qmin(n)*(a(n)-smin(n))

c	a(n)=a2(n)+a4(n)+a(n)/2.-wrig
	a(n)=a2(n)+a4(n)+a(n)/2.-wrig+
     *(b(n)-brig)*l(n)*dtdx

c ************ anti diffusion step ***********************
	do j=2,n-1
	dp(j)=ac(j+1)-ac(j)
	enddo

	do j=2,n-1
		if(dp(j)<0.)then
			dpl(j)=-1.
		elseif (dp(j)>0.)then 
			dpl(j)=1.
		else
			dpl(j)=0.
		endif
	enddo

	do j=2,n-1

	koef1=dpl(j)*dp(j-1)
	koef2=abs(dp(j))
	koef2=koef2/8.
	koef3=dpl(j)*dp(j+1)

	koef4=min1(koef1,koef2,koef3)
	koef5=max1(0.0,koef4)

	fc(j)=dpl(j)*koef5

	enddo

	fc(1)=0.
	fc(n)=0.
	fc(n-1)=0.
	
	do j=2,n-1
c	a(j)=ac(j)
	a(j)=ac(j)-(fc(j)-fc(j-1))
	enddo
c ************ end of anti diffusion step *********************** 
  
   	return
      end subroutine shasta
          
      end module utils